// Ring: pid 0 sends a counter to pid 1. Pid 1 increments it by pid.
// pid 1 sends to pid 2; pid2 increments by pid, and so on

// mpic++ ring.cpp
// mpirun -np 100 --oversubscribe. ./a.out

#include <iostream>
#include "mpi.h"
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <numeric>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  int pid = 0, nproc = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  int counter = 0;
  int tag = 0;
  
  //ring
  int src = (pid - 1 + nproc)%nproc;
  int dest = (pid + 1)%nproc; //%nproc lo vuelve periodico
  if (0 == pid){
    MPI_Send(&counter, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
    MPI_Recv(&counter, 1, MPI_INT, dest, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  } else{
    MPI_Recv(&counter, 1, MPI_INT, dest, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    counter += pid;
    MPI_Send(&counter, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
  }  

  if (0 == pid) {
    std::cout << counter << std::endl;
  }
  
  MPI_Finalize();
  return 0;
}

//MPI_Send(&counter, 1, MPI_INT, 1, tag, MPI_COMM_WORLD);
//MPI_Recv(&counter, 1, MPI_INT, 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);