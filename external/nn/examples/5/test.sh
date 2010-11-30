#!/bin/bash

N=373
M=455

if [ ! -x ../../nnbathy ]
then
    echo "error: no executable found"
    echo 'Run "./configure" and "make" in the source directory'
    exit 1
fi

echo ""
echo -n "Linear interpolation..."
../../nnbathy -i data.txt -n "$N"x"$M" -P alg=l > lin.txt
echo "done"
echo -n "Natural Neighbours Sibson interpolation..."
../../nnbathy -i data.txt -n "$N"x"$M" -W 0 > nn.txt
echo "done"
echo ""
echo 'To visualize, in Matlab run "viewexample"'
