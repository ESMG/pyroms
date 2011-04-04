#!/bin/bash

if [ ! -x ../../csabathy ]
then
    echo "error: no executable found"
    echo 'Run "./configure" and "make" in the upper level directory'
    exit 1
fi

echo "  Following are two examples of approximation of data with known error."

echo
echo "  The first data set consists of:"
echo "  - N random points of Franke test function with standard deviation of 0.2"
echo "  - N random points of linear function z = x + y with standard deviation of 5.0"
echo "  (N = 100, 300, 1000)."
echo

for N in 100 300 1000
do 
    echo "  N = $N:"
    echo -n "    Generating..."
    echo "1" | awk -f ./generate.awk N=$N F=0.2 P=5 > data-F-${N}.txt
    echo "done"
    echo -n "    Approximating into 256x256 grid..."
    ../../csabathy -i data-F-$N.txt -n 256x256 -P nppc=15 -P k=70 > out-F-${N}.txt
    echo "done"
done

echo
echo "  The second data set consists of:"
echo "  - N random points of Franke test function with standard deviation of 5.0"
echo "  - N random points of linear function z = x + y with standard deviation of 0.2"
echo "  (N = 100, 300, 1000)."
echo

for N in 100 300 1000
do 
    echo "  N = $N:"
    echo -n "    Generating..."
    echo "1" | awk -f ./generate.awk N=$N F=5 P=0.2 > data-P-${N}.txt
    echo "done"
    echo -n "    Approximating into 256x256 grid..."
    ../../csabathy -i data-P-$N.txt -n 256x256 -P nppc=15 -P k=70 > out-P-${N}.txt
    echo "done"
done

echo
echo "  Finished"
echo "  To visualize, in Matlab run \"viewexample\""
