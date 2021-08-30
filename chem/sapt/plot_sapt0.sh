#!/usr/bin/env bash

find . -mindepth 2 -name "*log" | 
xargs grep "Total sSAPT0" | 
sed 's/.\///;s/-sapt.*g://' | 
awk '{print $1,$(NF-1)}' | 
sort -k1n | 
gnuplot --persist -e "set ylabel 'SAPT0 [kJ/mol]'; plot '-' with lines notitle"
