set terminal eps
set output "neighbors.eps"

set ticslevel 0
set title "Searching of Neighbors"

set xlabel "X"
set ylabel "Y"
set zlabel "Z"

splot "NN_test_particle.output" u 2:3:4 w p lc rgb "red" pt 7 title "Particles", "NN_test_neighbors.output" u 2:3:4 w p lc rgb "blue" pt 7 ps 0.5 title "Neighbors"
