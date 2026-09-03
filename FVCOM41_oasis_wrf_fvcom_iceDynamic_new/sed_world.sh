#!/bin/bash

for INP in *.F
do
echo $INP
#sed -e s"/MPI_COMM_WORLD/il_commlocal/" $INP > a
sed -e s"/MPI_COMM_WORLD/il_commlocal/" $INP > a
#sed -e s"/il_commlocal/MPI_COMM_WORLD/" $INP > a
mv a $INP
done

for INP in *.F
do
echo $INP
#sed -e s"/mpi_comm_world/il_commlocal/" $INP > a
sed -e s"/mpi_comm_world/il_commlocal/" $INP > a
#sed -e s"/il_commlocal/mpi_comm_world/" $INP > a
mv a $INP
done
