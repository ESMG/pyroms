#!/bin/sh

LOCALDIR=/usr/local
CURDIR=`pwd`

echo
echo "installing pyroms..."
echo
python setup.py install
echo "installing external libraries..."
echo "installing gridgen..."
cd $LOCALDIR/python/pyroms/external/nn
./configure
make install
cd $LOCALDIR/python/pyroms/external/csa
./configure
make install
cd $LOCALDIR/python/pyroms/external/gridutils
CPPFLAGS=-I$LOCALDIR/include ./configure
make install
cd $LOCALDIR/python/pyroms/external/gridgen
CPPFLAGS=-I$LOCALDIR/include CFLAGS=-I$LOCALDIR/include ./configure
make
make lib
make shlib
make install
PYROMS_PATH=`python -c 'import pyroms ; print pyroms.__path__[0]'`
cp $LOCALDIR/lib/libgridgen.so $PYROMS_PATH
echo "installing scrip..."
cd $LOCALDIR/python/pyroms/external/scrip/source
awk '{gsub(/LIBDIR = \/opt\/local\/lib/,"LIBDIR = '$MACPORTDIR'/lib");print}' makefile > makefile2
awk '{gsub(/INCDIR = \/opt\/local\/include/,"INCDIR = '$MACPORTDIR'/include");print}' makefile2 > makefile
awk '{gsub(/PREFIX = \/usr\/local/,"PREFIX = '$LOCALDIR'");print}' makefile > makefile2
mv -f makefile2 makefile
make
make f2py
make install
cp $LOCALDIR/lib/scrip.so $PYROMS_PATH/remapping
cd $CURDIR
echo
echo "Done installing pyroms..."
echo "pyroms make use of the so-called gridid file to store"
echo "grid information like the path to the grid file, the"
echo "number of vertical level, the vertical transformation"
echo "use, ..."
echo "Please set the environment variable PYROMS_GRIDID_FILE"
echo "to point to your gridid file. A gridid file template"
echo "is available here:"
echo "$LOCALDIR/python/pyroms/pyroms/gridid.txt"
read -p "Press any key to continue or Ctrl+C to quit this install"
echo
