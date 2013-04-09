#
set terminal jpeg size 1100,400
set output 'plot-999.jpeg'

set multiplot

unset key 
set xtics nomirror
set ytics nomirror
set yr[0:0.1]
set size 1,1
set origin 0.00, 0.00
set xlabel("distance [A]")
set ylabel("Intensity [arb. u.]")
plot 'output.dat127' u 1:3 w l lw 2 


set yr[0:1]
unset key 
unset y2tics
unset xtics
#set xlabel ""
#set ylabel ""
set size 0.4,0.4
set origin  0.58, 0.58
plot 'output.dat999' u 1:2 w l lw 2 

exit multiplot
