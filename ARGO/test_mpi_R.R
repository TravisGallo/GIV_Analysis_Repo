library(Rmpi, lib.loc = "/home/hgallo/Rlibs/")

# detect the number of CPU's avaliable
no_cpus <- mpi.universe.size() - 1

print(no_cpus)

# set up cluster
cl <- makeCluster(no_cpus, type = "MPI")

print(cl)

## Tell all workers to close down
stopCluster(cl)