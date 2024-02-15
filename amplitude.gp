set logscale x
set grid
set title 'Amplitudengang'
set xlabel 'Frequenz (rad/s)'
set ylabel 'Amplitude (dB)'
set terminal png
set output 'amplitude.png'
unset key
plot 'bode_data_magnitude.txt' using 1:2 linecolor 1 with lines