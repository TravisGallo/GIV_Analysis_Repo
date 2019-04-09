

## ssh into ARGO
ssh hgallo@argo.orc.gmu.edu

## change directory to scratch
cd /scratch/hgallo

## get an allocation to test script
salloc --nodes=5 --ntasks=20 --mem=16GB

## load modules for Rmpi
module purge
module load R/3.5.1
module load gcc/7.3.1
module load openmpi/gcc/64/1.10.1
R

## had to get a build for Rmpi from ORC
install.packages("Rmpi", configure.args=c("--with-mpi=$MPI_HOME"))

## Install other R packages
install.packages(c("iterators","igraph", "sp", "raster", "doParallel", "gdistance", "crayon","dplyr", 
                   "data.table"))

install.packages("igraph")

