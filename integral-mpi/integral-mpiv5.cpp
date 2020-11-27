// mpic++ -g integral-mpiv5.cpp
// mpirun -np 8 ./a.out
// Wuilson Estacio
#include <cstdio>
#include <cstdlib>
#include "mpi.h"

/* Problem parameters */
double f(double x) {
  return x*x;
}

const int numberRects = 800000000;
const double lowerLimit = 2.0;
const double upperLimit = 5.0;

int main(int argc, char **argv) 
{
  int i; 
  double area;
  double at;
  double heigth;
  double width;

  int pid, np;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  int nrect = numberRects/np;
  int imin = pid*nrect;

  
  double tstart = MPI_Wtime();
  area = 0.0;
  width = (upperLimit - lowerLimit) / numberRects;

  for (i = imin; i < imin + nrect; ++i) {
    at = lowerLimit + i*width + width/2.0;
    heigth = f(at);
    area = area + width*heigth;
  }
  double tend = MPI_Wtime();
  std::printf("Time from pid %d: %lf\n", pid, tend-tstart);
  
  //std::printf("The area (pid %d) from %lf to %lf is : %le\n", pid, lowerLimit + imin*width, lowerLimit + (imin+nrect)*width, area);

  int tag = 0;
  if (0 == pid) {
    double total = area;
    for (int src = 1; src < np; ++src) {
      MPI_Recv(&area, 1, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      total += area;
    }
    std::printf("El area total es %lf\n", total);
  } else {
    int dest = 0;
    MPI_Send(&area, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
  }
  
  
  MPI_Finalize();
  
  return 0;
}