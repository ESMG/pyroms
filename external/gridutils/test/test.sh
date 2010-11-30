#!/bin/bash

if [ ! -x ../getnodes ]
then
    echo "error: no ../getnodes found"
    echo 'run "./configure" and "make" in the upper level directory'
    exit 1
fi

echo ""
echo "input double density grid: gridnodes.txt"
echo ""

echo "1. Validating and printing some stats:"
echo ""
../getnodes gridpoints_DD.txt -i DD -v > /dev/null
echo ""

echo -n "2. Extracting cell corner nodes from a double-density grid..."
../getnodes gridpoints_DD.txt -i DD -o CO > gridpoints_CO.txt
echo "done"
echo ""

echo -n "3. Getting the boundary and writing it to bound.txt..."
../getbound gridpoints_DD.txt > bound.txt
echo "done"
echo ""

echo -n "4. Getting the boundary in index space and writing it to bound-r.txt..."
../getbound gridpoints_DD.txt -r > bound-r.txt
echo "done"
echo ""

echo -n "5. As above, with compacting, writing to bound-c-r.txt..."
../getbound gridpoints_DD.txt -c -r > bound-c-r.txt
echo "done"
echo ""

echo "6. Converting a few points to and from index space:"
echo ""
echo -n '   513252.3881 5186890.274 -> '
echo "513252.3881 5186890.274" | ../xy2ij -g gridpoints_DD.txt -o stdin
echo "   and back:"
echo -n "   (index) "
echo "513252.3881 5186890.274" | ../xy2ij -g gridpoints_DD.txt -o stdin |tr -d "\n"
echo -n '-> '
echo "513252.3881 5186890.274" | ../xy2ij -g gridpoints_DD.txt -o stdin | ../xy2ij -g gridpoints_DD.txt -o stdin -r
echo ""
echo -n '   (index) 20.5 10.5 -> '
echo "20.5 10.5" | ../xy2ij -g gridpoints_DD.txt -o stdin -r
echo "   and back:"
echo -n "   "`echo "20.5 10.5" | ../xy2ij -g gridpoints_DD.txt -o stdin -r |tr -d "\n"`
echo -n '-> '
echo "20.5 10.5" | ../xy2ij -g gridpoints_DD.txt -o stdin -r | ../xy2ij -g gridpoints_DD.txt -o stdin
echo ""

../getnodes gridpoints_DD.txt -x | sed '1,1d' > x.txt
../getnodes gridpoints_DD.txt -y | sed '1,1d' > y.txt

if [ ! -x ../gridbathy ]
then
    echo "no ../gridbathy found"
    echo "omitting tests for gridbathy"
    exit 0
fi

echo -n "7. Interpolating bathymetry with bivariate cubic spline..."
../gridbathy -b bathy.txt -g gridpoints_DD.txt > bathy-cs.txt
echo "done"
echo ""

echo -n "8. Interpolating bathymetry with linear interpolation..."
../gridbathy -b bathy.txt -g gridpoints_DD.txt -a 3 > bathy-l.txt
echo "done"
echo ""

echo -n "9. Interpolating bathymetry with Natural Neighbours interpolation..."
../gridbathy -b bathy.txt -g gridpoints_DD.txt -a 2 > bathy-nn.txt
echo "done"
echo ""

echo -n "10. Interpolating bathymetry with Non-Sibsonian NN interpolation..."
../gridbathy -b bathy.txt -g gridpoints_DD.txt -a 1 > bathy-ns.txt
echo "done"
echo ""

echo "to visualise grid, in matlab run \"viewgrid('gridpoints_CO.txt')\""
echo "to visualise grid boundary, in matlab run \"viewbound('bound.txt')\""
echo "  (also \"viewbound('bound-r.txt')\", \"viewbound('bound-c-r.txt')\")"
echo "to visualise interpolated bathymetry, in matlab run \"viewbathy\""
echo "to visualise bathymetry data, in matlab run \"viewbathydata('bathy.txt', 60)\""




