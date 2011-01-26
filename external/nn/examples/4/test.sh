#!/bin/bash

N=201

if [ ! -x ../../nnbathy ]
then
    echo "error: no executable found"
    echo 'Run "./configure" and "make" in the source directory'
    exit 1
fi

echo ""
echo -n "Linear interpolation..."
../../nnbathy -i data.txt -n "$N"x"$N" -P alg=l > lin.txt
echo "done"
echo -n "Natural Neighbours Sibson interpolation..."
../../nnbathy -i data.txt -n "$N"x"$N" > nn-inf.txt
echo "done"
echo -n "Natural Neighbours Sibson interpolation with wmin = 0..."
../../nnbathy -i data.txt -n "$N"x"$N" -W 0 > nn-0.txt
echo "done"
echo -n "Natural Neighbours Non-Sibsonian interpolation with wmin = 0..."
../../nnbathy -i data.txt -n "$N"x"$N" -P alg=ns -W 0 > nn-ns.txt
echo "done"
echo ""
echo 'To visualize, in Matlab run "viewexample"'
