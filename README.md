# fast-indexer

A fast indexer based on [xgandalf](https://www.desy.de/~twhite/crystfel/manual-indexamajig.html).

To install, you first need to install ```ninja```, ```meson```, ```eigen3``` and a decent C++ compiler. On Arch Linux this can be done with:

```bash
sudo pacman -S ninja meson
```

Now let's compile it with: 

```bash
cd fast_indexer_dev
./compile-xgandalf.sh cci
./compile-test.sh cci
```

In order to compile you only need a decent C++ compilere and build tools installed. To test the code you can run:

```bash
./test.sh test/data/image0_peakfinder8.txt
```

In order to use it from anywhere and call it from python, lets first modify the bash script with the following command:

```bash
rdir=$(pwd)/ ; sed '0,/\${ROOT:-.}/s/\${ROOT:-.}/CHANGE/' test.sh | sed "s#CHANGE#$rdir#" > test.sh
```

Finally let's make the binary accessible globally linking to a proper global path:

```bash
fatstindexingdir=$(pwd)
cd ~/bin/
ln -s ${fatstindexingdir}/test.sh xgandalfHCS
```

## Example usage

To see an example of how to use it (on ra.psi.ch) see the [notebook](./example/Example.ipynb) in the folder example.
