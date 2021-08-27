set terminal pngcairo size 1280,1024 enhanced linewidth 2 font "Helvetica, 16"

set xrange [1.95 : 2]
do for [i=1:500] {
        set output sprintf("output/emissionLEATest/emissionManyOp_%06d.png", i)
        set multiplot layout 2,2 columnsfirst title sprintf("t = %.2f ns", (i-1)*0.1)
        # Electron and Ion density
        set ylabel "density m-3"
        set yrange [-1 : 5e5]
        #set log y
        set key bottom left
        plot sprintf("output/emissionLEATest/emissionManyOutputs_fluid_000%03.0f.txt", i)\
        using ($1*1e+3):($3/1e15) with lines lt 1 lw 3 title "LFA ne",\
        "" using ($1*1e+3):($4/1e15) with points ps 0.5 pt 4 title "LFA ni",\
        sprintf("output/emissionLEATest/emissionManyOutputsLEA_fluid_000%03.0f.txt", i) \
        using ($1*1e3):($3/1e15) with lines lt 1 lw 2 lc 6 title "LEA ne",\
        "" using ($1*1e3):($4/1e15) with points ps 0.5 pt 5 title "LEA ni",\
        sprintf("output/emissionLEATest/emission_particle_000%03.0f.txt", i) \
        using ($1*1e3):($3/1e15) with lines lt 1 lw 2 lc 2 title "part ne",\
        "" using ($1*1e3):($4/1e15) with points ps 0.5 pt 3 title "part ni"
        #unset yrange
        #unset log y
        # Electric  field
        set key top center
        set yrange [-1e6:4e7]
        set ylabel "E(V/m)"
        plot sprintf("output/emissionLEATest/emissionManyOutputs_fluid_000%03.0f.txt", i)\
        using ($1*1e+3):2 with lines lt 1 lw 3 title "LFA",\
        sprintf("output/emissionLEATest/emissionManyOutputsLEA_fluid_000%03.0f.txt", i) \
        using ($1*1e3):2 with lines lt 1 lw 2 lc 6 title "LEA",\
        sprintf("output/emissionLEATest/emission_particle_000%03.0f.txt", i) \
        using ($1*1e3):2 with lines lt 1 lw 2 lc 2 title "part"
        unset xlabel
        unset yrange
        # Net charge
        set ylabel "density m-3"
        set yrange [-1 : 1e5]
        #set log y
        set key bottom left
        plot sprintf("output/emissionLEATest/emissionManyOutputs_fluid_000%03.0f.txt", i)\
        using ($1*1e+3):($4-$3)/1e15 with lines lt 1 lw 3 title "LFA ne-ni",\
        sprintf("output/emissionLEATest/emissionManyOutputsLEA_fluid_000%03.0f.txt", i) \
        using ($1*1e3):($4-$3)/1e15 with lines lt 1 lw 2 lc 6 title "LEA ne-ni",\
        sprintf("output/emissionLEATest/emission_particle_000%03.0f.txt", i) \
        using ($1*1e3):($4-$3)/1e15 with lines lt 1 lw 2 lc 2 title "part ne-ni"
        unset yrange
        # Fluxes
        set ylabel "Fluxes"
        set yrange [-1e9 : 1e9]
        #set log y
        set key bottom left
        plot sprintf("output/emissionLEATest/emissionManyOutputs_fluid_000%03.0f.txt", i)\
        using ($1*1e+3):($9/1e15) with lines lt 1 lw 3  lc 1title "LFA advective",\
        "" using ($1*1e+3):($10/1e15) with lines lt 1 lw 3 lc 2 title "LFA diff",\
        sprintf("output/emissionLEATest/emissionManyOutputsLEA_fluid_000%03.0f.txt", i) \
        using ($1*1e3):($9/1e15) with lines lt 1 lw 2 lc 3 title "LEA advective",\
        "" using ($1*1e3):($10/1e15) with lines lt 1 lw 2 lc 4 title "LEA diff"
        unset yrange
        unset multiplot
}

