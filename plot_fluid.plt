# To run the comment we should do:
#    gnuplot -p "data.gp"
#set terminal postfile       (These commented lines would be used to )
#set output  "d2_plot.ps"    (generate a postscript file.            )

#set terminal gif anumate delay 100
#set ouput "foorbar.gif"

t = t + 0.5/100.0
outfile = sprintf('animation/fluid%.4d.png',i)

set terminal png
set output outfile
set ticslevel 0
fluido = sprintf("<sort -n -k 7 state_%.4d", i)
titulo = sprintf("%d",i)
#splot [-10:10][-10:10][-10:10] fluido u 2:3:4:10 w p palette ps 0.5 pt 7 t titulo
#splot [-10:10][-10:10][-10:10] fluido u 2:3:4:10 w p ps 0.5 pt 7 t titulo
#splot fluido u 2:3:4:10 w p ps 0.5 pt 7 t titulo
plot [-10:10][-10:10] fluido u 2:3:10 w p palette ps 0.5 pt 7 t titulo

#data = sprintf('state_%.4d.dat',i)

i = i + 1
if (t<end_time) reread;
