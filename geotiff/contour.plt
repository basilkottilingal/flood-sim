set contour base
unset surface
set view map

set cntrparam levels 10    # 10 contour levels

splot 'SampleGrid.dat' matrix with lines
