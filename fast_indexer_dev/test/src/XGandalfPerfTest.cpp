// Copyright (2019-2021) Paul Scherrer Institute
// GPL

#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <algorithm>
#include <cstdlib>
#include <sstream>
#include <chrono>
#include <memory>
#include "xgandalf/ExperimentSettings.h"
#include "xgandalf/Lattice.h"
#include "xgandalf/IndexerPlain.h"

namespace {
    // Input file format:
    // input = unit_cell coordinate_vector  # numbers are white space/line separated, comment lines can be anywhere
    // vec = x y z
    // unit_cell = vec(a) vec(b) vec(c)     # on the first line of the input file (not the reciprocal lattice here)
    // coordinate_vector = vec{20,40}       # on line per vec (min 20, max 40 vecs)
    // comment_line = # comment
    struct InputData final {
        static constexpr int max_spots = 20000;
        Eigen::Matrix3f unit_cell{};
        Eigen::Matrix3Xf spot_matrix{3, max_spots}; // will be resized later
    };

    struct OutputData final {
        std::vector<xgandalf::Lattice> lattices;
        std::vector<int> peaks;
    };

    // vector coordinates
    struct Coord final {
        float x;
        float y;
        float z;
    };

    // program arguments
    struct Args final {
        std::string prog_name;
        std::string input_file_name;
        int repetitions;
        bool fast;
    } args{};

    // print error message (nullptr for none) and exit
    [[noreturn]] void error(const char* msg)
    {
        if (msg != nullptr)
            std::cerr << "Error: " << msg << '\n';
        std::exit((EXIT_FAILURE));
    }

    // print usage and error message (nullptr for none) and exit
    [[noreturn]] void usage(const char* msg)
    {
        std::cerr << "Usage: " << args.prog_name << " <input_file_name> [repetitions=1] [fast]\n";
        error(msg);
    }

    // accepted options for help
    std::array<std::string,3> help_opt = {"help", "-h", "--help"};

    // get program arguments
    void parse_args(int argc, char *argv[])
    {
        args.prog_name = argv[0];
        if ((argc < 2) || (argc > 4))
            usage("Wrong number of args");
        if (std::find(cbegin(help_opt), cend(help_opt), std::string{argv[1]}) != cend(help_opt))
            usage(nullptr);
        args.input_file_name = argv[1];
        args.repetitions = 1;
        args.fast = false;
        bool rep_ok = false;
        for (int i=2; i<argc; i++) {
            if (! args.fast && (argv[i][0] == 'f')) {
                args.fast = true;
            } else if (! rep_ok) {
                std::istringstream reps(argv[i]);
                reps >> args.repetitions;
                if (! reps)
                    usage("Cannot parse repetitions");
                rep_ok = true;
            } else {
                usage("Duplicate repetition argument");
            }
        }
    }

    // remove white space at left/start
    void ltrim(std::string &s) {
        s.erase(begin(s), std::find_if(cbegin(s), cend(s), [](unsigned char ch) {
            return !std::isspace(ch);
        }));
    }

    // remove white space at right/end
    void rtrim(std::string &s) {
        s.erase(std::find_if(crbegin(s), crend(s), [](unsigned char ch) {
            return !std::isspace(ch);
        }).base(), end(s));
    }

    // remove leading and trailing white space
    void trim(std::string &s) {
        ltrim(s);
        rtrim(s);
    }

    // input stream for next line in file
    std::istringstream next_line(std::ifstream& fin)
    {
        for (std::string line; std::getline(fin, line);) {
            if (line[0] == '#')
                continue;
            trim(line);
            if (line.empty())
                continue;
            return std::istringstream{line};
        }
        return std::istringstream{};
    }

    // get vector coordinates from file
    Coord read_coordinates(std::istringstream& iss)
    {
        Coord c;
        if (! (iss >> c.x >> c.y >> c.z))
            error("can't read vec");
        return c;
    }

    // read input data from file
    InputData read_data(const std::string& input_file_name)
    {
        std::ifstream input_file(input_file_name);
        if (! input_file.is_open())
            usage("unable to read file");

        InputData in_data;

        {
            std::istringstream line = next_line(input_file);
            Coord a = read_coordinates(line);
            Coord b = read_coordinates(line);
            Coord c = read_coordinates(line);
            Eigen::Matrix3f A;
            A << a.x, b.x, c.x, a.y, b.y, c.y, a.z, b.z, c.z;
            in_data.unit_cell = A.inverse().transpose();
        }

        int n_spots;
        for (n_spots=0; true; n_spots++) {
            std::istringstream line = next_line(input_file);
            if (!line || line.str().empty())
                break;
            Coord spot = read_coordinates(line);
            if (n_spots >= in_data.spot_matrix.cols())
                error("too many spots");
            in_data.spot_matrix.col(n_spots) << spot.x, spot.y, spot.z;
        };
        in_data.spot_matrix.conservativeResize(3, n_spots);

        return in_data;
    }
}

int main (int argc, char *argv[])
{
    parse_args(argc, argv);

    std::cout << "Reading data ...\n";
    InputData in_data = read_data(args.input_file_name);
    xgandalf::ExperimentSettings settings(12400.f,
                                          0.1f,
                                          0.1f, 0.05f, 0.005f,
                                          in_data.unit_cell, 0.02f, 0.1f);
    xgandalf::IndexerPlain indexer(settings);

    if (args.fast) {
        indexer.setGradientDescentIterationsCount(xgandalf::IndexerPlain::GradientDescentIterationsCount::few);
        indexer.setSamplingPitch(xgandalf::IndexerPlain::SamplingPitch::standard);
        indexer.setMaxPeaksToUseForIndexing(30);
    }

    std::cout << "Calling xgandalf ...\n";
    OutputData out_data;
    std::chrono::high_resolution_clock clock;
    auto t1 = clock.now();
    for (int rep=0; rep<args.repetitions; rep++) {
        indexer.index(out_data.lattices, in_data.spot_matrix, out_data.peaks);
    }
    auto t2 = clock.now();

    std::cout << "Found " << out_data.lattices.size() << " lattices\n";
    for (int i=0; i<out_data.lattices.size(); i++) {
        std::cout << "Latice " << i << ":\n";
        std::cout << out_data.lattices[i] << '\n';
    }

    double iter_time = std::chrono::duration<double, std::milli>(t2 - t1).count();
    std::cout << "repetitions: " << args.repetitions << ", time: " << iter_time << ", " << iter_time / args.repetitions << " ms per repetition\n";
}
