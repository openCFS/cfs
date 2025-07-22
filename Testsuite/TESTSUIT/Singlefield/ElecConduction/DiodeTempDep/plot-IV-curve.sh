#!/bin/bash


gnuplot -e "
set term png;
set output 'u-i-curve.png';
set xlabel 'U / V';
set ylabel 'I / A';
plot  '< paste history/DiodeTempDep-ms2-elecPotential-node-8-anoPt.hist history/DiodeTempDep-ms2-elecPotential-node-12-catPt.hist history/DiodeTempDep-ms2-elecCurrent-surfRegion-cathode.hist ' 
using (\$4-\$2):(-\$6) w lp lw 2
title 'U-I curve';
"
