// mpic++ mpi-v1.cpp
// mpirun -np 5 --oversubscribe ./a.out

#include <iostream>
#include "mpi.h"
#include <cstdlib>
#include <cmath>
#include <chrono>

const int NS = 10;

void fill(double *data, int ns, int nslocal, int pid, int nproc);
void print(double *data, int ns, int nslocal, int pid, int nproc);

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  int pid = 0, nproc = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  const int NSlocal = NS/nproc;
  double * data = new double [NSlocal] {0.0};
  
  fill(data, NS, NSlocal, pid, nproc);
  print(data, NS, NSlocal, pid, nproc);

  delete [] data;
  MPI_Finalize();
  return 0;
}

void fill(double *data, int ns, int nslocal, int pid, int nproc)
{
  for (int ilocal = 0; ilocal < nslocal; ++ilocal) {
    data[ilocal] = 2*ilocal;
  }
}

void print(double *data, int ns, int nslocal, int pid, int nproc)
{
  std::cout << "pid: " << pid << "\n";
  for (int ilocal = 0; ilocal < nslocal; ++ilocal) {
    std::cout << data[ilocal] << "  ";
  }
  std::cout << "\n";
}