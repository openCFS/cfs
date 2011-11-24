set datafile separator ","
set logscale y
set xlabel 'Frequency [Hz]'
set ylabel 'Lighthill Source (Amplitude) [??]'
plot 'acourhs_harm_p1.csv' using ($7):($3) with lines lw 2 title 'P1', \
     'acourhs_harm_p2.csv' using ($7):($3) with lines lw 2 title 'P2',
