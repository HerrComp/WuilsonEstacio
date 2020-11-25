// mpic++ -g integral-mpiv4.cpp
// mpirun -np 8 ./a.out
// Wuilson Estacio
#include <cstdio>
#include <cstdlib>
#include "mpi.h"

/* Problem parameters */
double f(double x) {
  return x*x;
}

const double lowerLimit = 2.0;
const double upperLimit = 6.0;

int main(int argc, char **argv) 
{
  std::cout.precision(15); 
  std::cout.setf(std::ios::scientific);

  int i; 
  double area;
  double at;
  double heigth;
  double width;

  int pid, np;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  int numberRects = 0;
  if (0 == pid) {
    std::printf("Escriba el numero de rectangulos en los que desea dividir la integral:\n");
    std::scanf("%d", &numberRects);
    std::printf("%d\n", numberRects);
  }
  MPI_Bcast(&numberRects, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
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
  
  std::printf("The area (pid %d) from %lf to %lf is : %le\n", pid, lowerLimit + imin*width, lowerLimit + (imin+nrect)*width, area);

  double total = 0;
  MPI_Reduce(&area, &total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  if (0 == pid) {
    std::printf("El area total de la int es %lf\n", total);
  }
  
  
  MPI_Finalize();
  
  return 0;
}