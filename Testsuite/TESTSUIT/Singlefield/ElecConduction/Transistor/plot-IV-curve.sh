#!/bin/bash


gnuplot -e "
set term png;
set output 'u-i-curve.png';
set xlabel 'Vgs / V';
set ylabel 'I / A';
plot  '< paste history/Transistor-ms1-elecPotential-node-12-drainPt.hist history/Transistor-ms1-elecPotential-node-8-sourcePt.hist history/Transistor-ms1-elecPotential-node-23-gatePt.hist history/Transistor-ms1-elecCurrent-surfRegion-Udrain.hist' 
using (\$6-\$4):(-\$8) w lp lw 2
title 'U-I curve';
"
