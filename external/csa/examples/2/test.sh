#!/bin/bash

if [ ! -x ../../csabathy ]
then
    echo "error: no executable found"
    echo 'Run "./configure" and "make" in the upper level directory'
    exit 1
fi

echo -n "  Approximating on 256x256 grid using default settings..."
../../csabathy -i data.txt -n 256x256 > out_default.txt
echo "done"
echo -n "  Approximating on 256x256 grid using -P nppc=20..."
../../csabathy -i data.txt -n 256x256 -P nppc=20 > out_nppc20.txt
echo "done"
echo -n "  Approximating on 256x256 grid using -P k=70..."
../../csabathy -i data.txt -n 256x256 -P k=70 > out_k70.txt
echo "done"
echo -n "  Approximating on 256x256 grid using -P k=280..."
../../csabathy -i data.txt -n 256x256 -P k=280 > out_k280.txt
echo "done"
echo "  Finished"
echo "  To visualize, in Matlab run \"viewexample\""
