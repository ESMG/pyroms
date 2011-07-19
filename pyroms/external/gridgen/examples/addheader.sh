#!/bin/bash

sed -n '$p' $1 > header.tmp
cat header.tmp $2 > "$2".tmp
mv -f "$2".tmp "$2"
rm -f header.tmp
