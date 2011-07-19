#!/bin/bash

if [ ! -x ../../nnbathy ]
then
    echo "error: no executable found"
    echo 'Run "./configure" and "make" in the source directory'
    exit 1
fi

echo "  Example of Natural Neighbours interpolation:"
echo "  Franke test function reconstruction by 100, 300 and 1000 random points"

for N in 100 300 1000
do 
    echo "  N = $N:"
    echo -n "    Generating..."
    echo "1" | awk -f ./generate.awk N=$N > data-${N}.txt
    echo "done"
    echo -n "    Interpolating into 256x256 grid..."
    ../../nnbathy -i data-$N.txt -n 256x256 > out-${N}.txt
    echo "done"
done

echo "  Finished"
echo "  To visualize, in Matlab run \"viewexample\""
