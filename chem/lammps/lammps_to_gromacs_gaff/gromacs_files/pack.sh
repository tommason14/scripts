#!/usr/bin/env bash

box=12
# 20 chains with 24 sulfonates each
gmx_mpi insert-molecules -ci polymer.gro -nmol 20 -box $box $box $box -try 200 -o pack.gro
# 10 waters per sulfonate
gmx_mpi insert-molecules -f pack.gro -ci tip4p.gro -nmol 4800 -box $box $box $box -try 200 -o pack.gro
# 0.2 M NaCl 
gmx_mpi insert-molecules -f pack.gro -ci cl.gro -nmol 96 -box $box $box $box -try 200 -o pack.gro
# 480 sodiums to balance polymer charge, and 96 for 0.2 M NaCl
gmx_mpi insert-molecules -f pack.gro -ci sodium.gro -nmol 576 -box $box $box $box -try 200 -o pack.gro
