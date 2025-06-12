set yrange [0:0.9909636]
set arrow 1  from 5,0 to 5,0.825803 nohead lt rgb "gray" lw 2
set arrow 2  from 10,0 to 10,0.825803 nohead lt rgb "gray" lw 2
set arrow 3  from 15,0 to 15,0.825803 nohead lt rgb "gray" lw 2
unset ylabel
unset xlabel
set style rectangle back fc rgb "gray"
plot"BlochMode3D_sorted.bloch.dat" u 5 t "1. mode"  with linespoints  ,\
    "BlochMode3D_sorted.bloch.dat" u 6 t "2. mode"  with linespoints  ,\
    "BlochMode3D_sorted.bloch.dat" u 7 t "3. mode"  with linespoints  ,\
    "BlochMode3D_sorted.bloch.dat" u 8 t "4. mode"  with linespoints  ,\
    "BlochMode3D_sorted.bloch.dat" u 9 t "5. mode"  with linespoints 
replot
set object 1 rect from 5,0.450291 to 10,0.544357
set object 2 rect from 10,0.536999 to 15,0.615464
set object 3 rect from 0,0.544357 to 5,0.552298
set object 4 rect from 0,0.583002 to 5,0.586142
set object 5 rect from 10,0.702702 to 15,0.709994
