#!/bin/bash

if [ ! -x ../../nnbathy ]
then
    echo "error: no executable found"
    echo 'Run "./configure" and "make" in the source directory'
    exit 1
fi

echo ""
echo -n "Natural Neighbours Sibson interpolation with -W 0 ..."
../../nnbathy -W 0 -n 152x114 -x 591020.57127923 591321.93673645 -y 4260093.61151167 4259867.85217794 -i data.txt > nn.txt
echo "done"
echo ""
echo 'To visualize, in Matlab run "viewexample"'
