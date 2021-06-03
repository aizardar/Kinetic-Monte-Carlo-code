     unset surface
     set contour
     set cntrparam ...
     set table 'filename'
     splot ...
     unset table
  
     set term <whatever>
     plot 'filename'





set dgrid3d 50,50 gauss 10e-6,10e-6
set logscale z
splot "fort.12" u 1:2:7 w l
