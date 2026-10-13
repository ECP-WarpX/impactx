set xtics font  'helvetica,25'
set ytics font 'helvetica,25'
set xlabel 'x/\sigmax' font 'helvetica,30' offset 0,-1
set ylabel 'Ex (MV/m)' font 'helvetica,30' offset -1,0
set lmargin 12
set bmargin 6
set xrange [-6:6]
set yrange [-10:10]
set nokey
set mxtics 5
set mytics 5
f(x)=8.0*exp(-x**2/2.0)
plot f(x) w filledcu y1=0.0 lt 3 fs transparent solid 0.25
replot 'Ex_Mayes.dat' u 1:2 w l lw 2 lt 1
replot 'Ex_Mayes.dat' u 1:3 w l lw 2 lt 6
replot 'Ex_Mayes.dat' u 1:4 w l lw 2 lt 8
replot 'Ex_Mayes.dat' u 1:5 w l lw 2 lt 7
