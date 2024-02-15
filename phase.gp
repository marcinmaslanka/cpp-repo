set logscale x
set grid
set title 'Phasengang'
set xlabel 'Frequenz (rad/s)'
set ylabel 'Phase (grad)'
set terminal png
set output 'phase.png'
unset key
plot 'bode_data_phase.txt' using 1:2 linecolor 1 with lines