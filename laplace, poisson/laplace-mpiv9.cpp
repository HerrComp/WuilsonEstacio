/*
//modela arreglo unidim
typedef std::vector<double> Matrix;

// mpic++ -std=c++17 -fsanitize=address -fconcepts -g -O3 laplace-mpivf.cpp
//-Werror -Wall
// mpic++ -std=c++17 -O3 laplace-mpiv8.cpp
// mpirun -np 4 ./a.out 12
//Wuilson Estacio
/* gnuplot
set pm3d; set contour base
set term pdf; set out 'matrix.pdf'
splot 'datos1.txt' w pm3d
*/

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include "mpi.h"

int N = 10;
const double L = 1.479;
double DELTA = L/N;  
const int STEPS = 200;



// MPI 
void print_matrix_slice(double * array, int nx, int ny);
void mpi_print_matrix(int pid, int np, double * array, int nx, int ny); 
void mpi_interchange_data(int pid, int np, double * array, int nx, int ny); // nx = Nl + 2
void mpi_initial_conditions(int pid, int np, double * array, int nx, int ny); 
void mpi_boundary_conditions(int pid, int np, double * array, int nx, int ny); 
void mpi_evolve(int pid, int np, double * array, int nx, int ny); // 
void mpi_print_gnuplot(int pid, int np, double * array, int nx, int ny);
void mpi_print_gnuplot_slice(int pid, int np, double * array, int nx, int ny); 

int main(int argc, char **argv) {
  
  N = std::atoi(argv[1]);
  DELTA = L/N; 

  
  int pid, np;
  MPI_Init(&argc, &argv); 
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  
  int Nl = N/np; 
  double * data = new double [(Nl+2)*N] {0.0};
  
  std::fill(data, data + (Nl+2)*N, pid); 
  mpi_initial_conditions(pid, np, data, Nl+2, N);
  mpi_interchange_data(pid, np, data, Nl+2, N);
  mpi_boundary_conditions(pid, np, data, Nl+2, N);

  mpi_interchange_data(pid, np, data, Nl+2, N);

  for (int istep = 0; istep < STEPS; ++istep) {  
    mpi_evolve(pid, np, data, Nl+2, N);
    mpi_interchange_data(pid, np, data, Nl+2, N);
  }
  // imprimir
  mpi_print_matrix(pid, np, data, Nl + 2, N); 
  mpi_print_gnuplot(pid, np, data, Nl+2, N);
 
  delete [] data; 
  MPI_Finalize();
  return 0;
}

void mpi_print_matrix(int pid, int np, double * array, int nx, int ny){
  if (0 == pid) {
    print_matrix_slice(array, nx, ny);
    double * buffer = new double [nx*ny];
    for (int ipid = 1; ipid < np; ++ipid) {
      MPI_Recv(buffer, nx*ny, MPI_DOUBLE, ipid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      print_matrix_slice(buffer, nx, ny);
    }
    delete [] buffer; 
  } else {
    MPI_Send(array, nx*ny, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
}


void print_matrix_slice(double * array, int nx, int ny) { 
  for (int ii = 0 ; ii < nx; ++ii) {
    if (0 == ii or ii == nx-1) { 
      std::cout << "G : " ; 
    } else {
      std::cout << "  : " ;
    } 
    for (int jj = 0 ; jj < ny; ++jj) {
      std::printf("%6.1lf ", array[ii*ny + jj]);
    }  
    std::cout << "\n";
  }
}

void mpi_interchange_data(int pid, int np, double * array, int nx, int ny) { // nx = Nl + 2
  if (0 != pid) { 
    MPI_Send(array + ny, ny, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD);
    MPI_Recv(array, ny, MPI_DOUBLE, pid-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  if (np-1 != pid) {  
    MPI_Recv(array + ny*(nx - 1), ny, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(array + ny*(nx - 2), ny, MPI_DOUBLE, pid+1, 1, MPI_COMM_WORLD);
  }
}

void mpi_initial_conditions(int pid, int np, double * array, int nx, int ny) {
  for (int ii = 1; ii < nx-1; ++ii) {
    for (int jj = 0; jj < ny; ++jj) {
      array[ii * ny + jj] = 1.0;
    }
  }
}

/* 
(x,y), (x,L)=100, (L,y)=0, (0,y)=100, (x,0)=0,
*/
void mpi_boundary_conditions(int pid, int np, double * array, int nx, int ny) {
  int ii = 0, jj = 0;
  // (x,L)
  if (0 == pid) { 
    ii = 1; 
    for (jj = 0; jj < ny; ++jj)
      array[ii * ny + jj] = 100;
  }
  // (x,0)
  if (pid == np-1) { 
    ii = nx - 2; 
    for (jj = 0; jj < ny; ++jj)
      array[ii * ny + jj] = 0;
  }
  // (0,y)  
  jj = 0;
  for (ii = 1; ii < nx - 1; ++ii)
    array[ii * ny + jj] = 100;
  // (L,y)
  jj = ny - 1;
  for (ii = 1; ii < nx - 1; ++ii)
    array[ii * ny + jj] = 0;
}

 
void mpi_evolve(int pid, int np, double * array, int nx, int ny) {
  for (int ii = 1; ii < nx-1; ++ii) { 
    for (int jj = 0; jj < ny; ++jj) { 
      if (0 == pid && ii == 1) 
        continue;
      if (pid == np-1 && ii == nx - 2)  
        continue;
      if (jj == 0) 
        continue;
      if (jj == ny - 1) 
        continue;
      
      array[ii * ny + jj] = (array[(ii + 1) * ny + jj] + array[(ii - 1) * ny + jj] +
			     array[ii * ny + jj + 1] + array[ii * ny + jj - 1]) / 4.0;
    }
  }
}

void mpi_print_gnuplot(int pid, int np, double * array, int nx, int ny)
{
  if (0 == pid) {
    mpi_print_gnuplot_slice(pid, np, array, nx, ny); 
    double * buffer = new double [nx*ny]; 
    for (int ipid = 1; ipid < np; ++ipid) {
      MPI_Recv(buffer, nx*ny, MPI_DOUBLE, ipid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
      mpi_print_gnuplot_slice(ipid, np, buffer, nx, ny); 
    }
    delete [] buffer; 
  } else {
    MPI_Send(array, nx*ny, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
}

void mpi_print_gnuplot_slice(int pid, int np, double * array, int nx, int ny) {
  for (int ii = 1; ii < nx-1; ++ii) { 
    for (int jj = 0; jj < ny; ++jj) { 
      std::cout << (ii-1 + pid*(nx-2)) * DELTA << " " << jj * DELTA << " " << array[ii * ny + jj] << "\n"; 
    }
    std::cout << "\n"; 
  }
}



