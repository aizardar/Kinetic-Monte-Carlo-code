################################################
#			PS(greyscale)
################################################

set terminal postscript enhanced color 
set output 'plot_gs_600K.ps'


set view map
unset surface
set pm3d 
set palette color
set size ratio 2
set contour
set key outside
set cntrparam cubicspline

set cntrparam levels 5
set xrange[0:10]
set yrange[0:30]
set xlabel "Al Concentration (at.%)"
set ylabel "Mn Concentration (at.%)"

set dgrid3d 100,100 gauss 1

set palette grey

 


splot 'D_vs_Al_concentration_250K.dat' using 1:2:(log10($3)) title "Diffusivity(m^2/s)" with lines lt 2 


################################################################
#				PNG(greyscale)
################################################################
set terminal png enhanced size 1024,768
set output 'plot_gs_600K.png'


set view map
unset surface
set pm3d 
set palette color
set size ratio 2
set contour
set key outside
set cntrparam cubicspline

set cntrparam levels 5
set xrange[0:10]
set yrange[0:30]
set xlabel "Al Concentration (at.%)"
set ylabel "Mn Concentration (at.%)"

set dgrid3d 100,100 gauss 1

set palette grey


splot 'D_vs_Al_concentration_250K.dat' using 1:2:(log10($3)) title "Diffusivity(m^2/s)" with lines lt 2 


#################################################################
#			PNG(color)
#################################################################
set terminal png enhanced size 1024,768
set output 'plot_color_600K.png'


set view map
unset surface
set pm3d 
set palette color
set size ratio 2
set contour
set key outside
set cntrparam cubicspline

set cntrparam levels 5
set xrange[0:10]
set yrange[0:30]
set xlabel "Al Concentration (at.%)"
set ylabel "Mn Concentration (at.%)"

set dgrid3d 100,100 gauss 1

set palette model RGB
set palette defined
splot 'D_vs_Al_concentration_250K.dat' using 1:2:(log10($3)) title "Diffusivity(m^2/s)" with lines lt 2 

###################################################################
#			PS(color)
###################################################################
set terminal postscript enhanced color
set output 'plot_color_600K.ps'


set view map
unset surface
set pm3d 
set palette color
set size ratio 2
set contour
set key outside
set cntrparam cubicspline

set cntrparam levels 5
set xrange[0:10]
set yrange[0:30]
set xlabel "Al Concentration (at.%)"
set ylabel "Mn Concentration (at.%)"

set dgrid3d 100,100 gauss 1

set palette model RGB
set palette defined
splot 'D_vs_Al_concentration_250K.dat' using 1:2:(log10($3)) title "Diffusivity(m^2/s)" with lines lt 2 
