reset;
set encoding iso_8859_1
cd 'Ihr Verzeichnis';
set terminal x11;
pi = 4*atan(1)

# Initialize max variable with maximum of column vector
real_max=`cat directivity.dat | cut -d' ' -f7 | sed 's/[eE]/*10^/' | sed 's/+//' | sed 's/^/scale=12; /' | bc | sort -n | tail -1`
max=75
unit='dB'

set nokey
set size square
set parametric
unset xlabel
unset ylabel
set title "Polardarstellung" font "Times,24"
unset logscale x
unset logscale y
set xrange[-1.1:1.1]
set yrange[-1.1:1.1]
set samples 1500

set label 1 "0°" at 1.05*cos(0*pi)-0.025,1.05*sin(0*pi) font "Times,10"
set label 2 "45°" at 1.05*cos(0.25*pi)-0.025,1.05*sin(0.25*pi) font "Times,10"
set label 3 "90°" at 1.05*cos(0.5*pi)-0.025,1.05*sin(0.5*pi) font "Times,10"
set label 4 "135°" at 1.05*cos(0.75*pi)-0.025,1.05*sin(0.75*pi) font "Times,10"
set label 5 "180°" at 1.05*cos(1*pi)-0.025,1.05*sin(1*pi) font "Times,10"
set label 6 "-135°" at 1.05*cos(1.25*pi)-0.025,1.05*sin(1.25*pi) font "Times,10"
set label 7 "-90°" at 1.05*cos(1.5*pi)-0.025,1.05*sin(1.5*pi) font "Times,10"
set label 8 "-45°" at 1.05*cos(1.75*pi)-0.025,1.05*sin(1.75*pi) font "Times,10"

set label 10 sprintf("%.2f %s", 0.0*max, unit) at 1.05*0*cos(0.4*pi),1.05*0*sin(0.4*pi) font "Times 10"
set label 11 sprintf("%.2f %s", 0.2*max, unit) at 1.05*0.2*cos(0.4*pi),1.05*0.2*sin(0.4*pi) font "Times,10"
set label 12 sprintf("%.2f %s", 0.4*max, unit) at 1.05*0.4*cos(0.4*pi),1.05*0.4*sin(0.4*pi) font "Times,10"
set label 13 sprintf("%.2f %s", 0.6*max, unit) at 1.05*0.6*cos(0.4*pi),1.05*0.6*sin(0.4*pi) font "Times,10"
set label 14 sprintf("%.2f %s", 0.8*max, unit) at 1.05*0.8*cos(0.4*pi),1.05*0.8*sin(0.4*pi) font "Times,10"
set label 15 sprintf("%.2f %s", 1.0*max, unit) at 1.05*1*cos(0.4*pi),1.05*1*sin(0.4*pi) font "Times,10"

set label 3000 "Angle (deg)" at -0.8,0.8 rotate by 45 font "Times,12"
# set label 3001 "Reflexionsfaktor r" at 0.4*cos(0.42*pi),0.4*sin(0.42*pi) font "Times,12" rotate by 0.45*180

a=-1
set multiplot
load 'polarloop.gnuplot'
set trange [-pi:pi]
plot 'directivity.dat' using ($6/max*cos($3)):($6/max*sin($3)) with lines lw 3 title 'conforming'
unset multiplot
unset parametric
