reset
# set term wxt 501
set xlabel 'eps33, \%'
set ylabel 'Omega_11'
# set size square
# set yrange [] reverse
# set logscale x
# set format x "%g"
p 'result.out' u (-$4*100):($8) w l 
