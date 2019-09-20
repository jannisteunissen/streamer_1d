#!/bin/bash

if (( $# != 1 )); then
    echo "Usage example: $0 output/sim_fluid"
    exit 1
fi

mkdir -p pngs

bname="$1"
sim_name="${1##*/}"
png_name="pngs/$sim_name"
sfile="${bname}_scalars.txt"

lbound_elec="1e15"
ubound_field=$(awk 'NR==2 {max = $3 + 0; next} {if ($3 > max) max = $3;} END {printf "%.2e", max*1.1}' "$sfile") || exit 1
ubound_elec=$(awk 'NR==2 {max = $9 + 0; next} {if ($9 > max) max = $9;} END {printf "%.2e", max*2}' "$sfile") || exit 1
dt_output=$(awk 'NR==2 {t0 = $1 + 0}; NR==3 {print $1 - t0; exit 0}' "$sfile") || exit 1
n_files=$(ls "${bname}_"??????.txt | wc -l) || exit 1

gnuplot << EOF
set terminal pngcairo size 800,600 enhanced linewidth 2 font "Helvetica,16"
set key samplen 2
set key autotitle columnhead
set style data lines
dt = $dt_output

do for [i=1:$n_files] {
set key bottom center
set output sprintf("${png_name}_%06d.png", i)
set multiplot layout 2,1 title sprintf('t = %.2f ns', (i-1)*dt*1e9)
set ylabel "density (m^{-3})"
set yrange [$lbound_elec:$ubound_elec]
set logscale y
plot sprintf("${bname}_%06d.txt", i) u (\$1*1e3):3 title "n_e", "" u (\$1*1e3):4 title "n_i";
unset logscale y
set key top center
set ylabel "E (V/m)"
set yrange [-0.08*$ubound_field:$ubound_field]
set xlabel "x (mm)" offset 0,0.8
plot sprintf("${bname}_%06d.txt", i) u (\$1*1e3):2 title "E"
unset xlabel
unset multiplot
}
EOF

ffmpeg -y -i ${png_name}_%06d.png "${sim_name}.mp4"
