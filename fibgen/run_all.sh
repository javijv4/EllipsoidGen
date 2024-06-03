#!/bin/bash

for meshdir in '../mesh_coarse' '../mesh_mid' '../mesh_fine'
do
    cheartsolver.out laplace_trans.P -\#meshdir=${meshdir} --prep
    mpirun -np 5 cheartsolver.out laplace_trans.P  -\#meshdir=${meshdir}
    mpirun -np 5 cheartsolver.out laplace_long.P  -\#meshdir=${meshdir}
done