w1=0.95
w2=0.05
bg=0.0146003
Min=0.0
Max=0.025
Np=100000001
n=300
width=(Max-Min)/n
bin(x,width)=width*(floor((x-Min)/width)+0.5)+Min
r(x,y,z)=sqrt(x**2+y**2+z**2)
set table 'InitialBeam.dat'
plot 'beam_initial.csv' u (bin(r(bg*$6,$7,$8),width)):(1.0/(Np*width)) smooth freq w l lt 7
unset table
set table 'FinalBeam.dat'
plot 'beam_final.csv' u (bin(r(bg*$6,$7,$8),width)):(1.0/(Np*width)) smooth freq w l lt 7
unset table
set logscale y
set xrange [Min:Max]
set yrange [0.1:1.6e6]
set xtics font 'helvetica,25'
set ytics font 'helvetica,25'
set lmargin 18
set rmargin 8
set bmargin 5
set xlabel 'r (m)' font 'helvetica,30' offset 0,-1
set pointsize 0.7
set key font 'helvetica,20'
set grid
show grid
plot 'Pdensity.out.0' u 1:2 w l lw 2 lt 7 title 'Ideal density'
replot 'InitialBeam.dat' u 1:($2/(4.0*pi*($1)**2)) w l lw 2 lt 6 title 'Initial beam'
replot 'FinalBeam.dat' every 5 u 1:($2/(4.0*pi*($1)**2)) w p lt 8 title 'Final beam'
