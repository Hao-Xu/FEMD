reset
# set term wxt 501
set xlabel 'eps33, \%'
set ylabel 'stress33-stress11, MPa'
# set size square
# set yrange [] reverse
# set logscale x
# set format x "%g"
p 'result.out' u (-$4*100):(-($7-$5)/1e6) w l title 'FLAC3D',\
'../UMAT/result.out' u (-$4*100):(-($7-$5)/1e6) w l t 'ABAQUS',\
'../UMAT2/result.out' u (-$4*100):(-($7-$5)/1e6) w l t 'DYNAFLOW'
