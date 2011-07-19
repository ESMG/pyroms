#!/bin/bash

if [ ! -x ../../nnbathy ]
then
    echo "error: no executable found"
    echo 'Run "./configure" and "make" in the source directory'
    exit 1
fi

echo | awk -f generate-data.awk > data.txt
echo | awk -f generate-points.awk > points.txt

../../nnbathy -i data.txt -o points.txt > results.txt
../../nnbathy -i data.txt -o points.txt -P alg=ns > results-ns.txt

paste -d" " results.txt results-ns.txt points.txt |\
cut -f"1 2 3 6 9" -d" " |\
awk '{a = $3 - $5; s1 += (a > 0) ? a : -a; b = $4 - $5; s2 += (b > 0) ? b : -b;} END {print "  sum(|z_sibson_i - z_i|) =", s1; print "  sum(|z_nonsibson_i - z_i|) =", s2;}'
