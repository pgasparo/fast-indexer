# fast-indexer

A fast indexer based on [xgandalf](https://www.desy.de/~twhite/crystfel/manual-indexamajig.html).

### Install

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

#### Quick fix in case installation fails... 

In case ```compile-test.sh``` is raising errors, go to the back to the ```fast_indexer_dev``` folder folder and create a new directory:

```bash
cd fast_indexer_dev
mkdir test-build
```

Then get the proper paths:

```bash
cd test/src
PATHI=$(pkg-config --cflags /tmp/fast-indexer/fast_indexer_dev/xgandalf-install/lib/pkgconfig/xgandalf.pc)
# OUTPUT: 
# -I/tmp/fast-indexer/fast_indexer_dev/xgandalf-install/include -DEIGEN_NO_DEBUG -DEIGEN_NO_AUTOMATIC_RESIZING -I/usr/include/eigen3 
PATHL=$(pkg-config --libs /tmp/fast-indexer/fast_indexer_dev/xgandalf-install/lib/pkgconfig/xgandalf.pc)
# OUTPUT:
# -L/tmp/fast-indexer/fast_indexer_dev/xgandalf-install/lib -lxgandalf
```

now from you ```test/src/``` folder compile with the command:

```bash
g++ $PATHI -DEIGEN_NO_DEBUG -DEIGEN_NO_AUTOMATIC_RESIZING -I/usr/include/eigen3 -o ../../test-build/xgandalf XGandalfPerfTest.cpp $PATHL
```

Finally, we need to set correctly the library path for xgandal:

```bash
export LD_LIBRARY_PATH=/tmp/fast-indexer/fast_indexer_dev/xgandalf-install/lib/
```

Don't forge to update the path with yours.

### Test

To test the code you can run:

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
