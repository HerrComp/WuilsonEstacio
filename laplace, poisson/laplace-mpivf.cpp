// mpic++ -std=c++17 -fsanitize=address -fconcepts -g -o3 laplace-mpivf.cpp
//-Werror -Wall
// mpirun -np 4 ./a.out
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include "mpi.h"

// constants
//const int N = int(L / DELTA) + 1;
int N = 10;
const double L = 1.479;
double DELTA = L/N;
const int STEPS = 200;

typedef std::vector<double> Matrix;

void initial_conditions(Matrix &m);
void boundary_conditions(Matrix &m);
void evolve(Matrix &m);
void print_gnuplot(const Matrix &m);
void print_matrix(const Matrix &m);
void init_gnuplot(void);
void plot_gnuplot(const Matrix &m);

////////////////// MPI versions
void print_matrix_slice(double * array, int nx, int ny);
void mpi_print_matrix(int pid, int np, double * array, int nx, int ny);
void mpi_interchange_data(int pid, int np, double * array, int nx, int ny); // nx = Nl + 
void mpi_initial_conditions(int pid, int np, double * array, int nx, int ny);
void mpi_boundary_conditions(int pid, int np, double * array, int nx, int ny);
void mpi_evolve(int pid, int np, double * array, int nx, int ny);
void mpi_print_gnuplot(int pid, int np, double * array, int nx, int ny);
void mpi_print_gnuplot_slice(int pid, int np, double * array, int nx, int ny);

int main(int argc, char **argv) {
  /*
  Matrix data(N * N);
  initial_conditions(data);
  boundary_conditions(data);

  //init_gnuplot();
  for (int istep = 0; istep < STEPS; ++istep) {
    evolve(data);
    //plot_gnuplot(data);
  }
  //print_gnuplot(data);
  print_matrix(data);
  */
  
  N = std::atoi(argv[1]);
  DELTA = L/N;

  int pid, np;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  // local array
  int Nl = N/np;
  double * data = new double [(Nl+2)*N] {0.0};
  // llenar con el pid
  std::fill(data, data + (Nl+2)*N, pid); 
  // initial and boundary conditions
  mpi_initial_conditions(pid, np, data, Nl+2, N);
  mpi_interchange_data(pid, np, data, Nl+2, N);
  mpi_boundary_conditions(pid, np, data, Nl+2, N);
  mpi_interchange_data(pid, np, data, Nl+2, N);
  // evolution
  for (int istep = 0; istep < STEPS; ++istep) {  
    mpi_evolve(pid, np, data, Nl+2, N);
    mpi_interchange_data(pid, np, data, Nl+2, N);
  }
  // imprimir
  //mpi_print_matrix(pid, np, data, Nl + 2, N);
  //mpi_print_gnuplot(pid, np, data, Nl+2, N);
  
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
      //std::cout << array[ii*ny + jj] << "  ";
      std::printf("%5.1lf ", array[ii*ny + jj]);
    }
    std::cout << "\n";
  }
}

void mpi_interchange_data(int pid, int np, double * array, int nx, int ny) { // nx = Nl + 2
  if (0 != pid) {
    //enviar hacia atras
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

void mpi_boundary_conditions(int pid, int np, double * array, int nx, int ny) {
  int ii = 0, jj = 0;

  if (0 == pid) {
    ii = 1; // fila real
    for (jj = 0; jj < ny; ++jj)
      array[ii * ny + jj] = 100;
  }

  if (pid == np-1) {
    ii = nx - 2; // fila real
    for (jj = 0; jj < ny; ++jj)
      array[ii * ny + jj] = 0;
  }
  
  jj = 0;
  for (ii = 1; ii < nx - 1; ++ii)
    array[ii * ny + jj] = 0;

  jj = ny - 1;
  for (ii = 1; ii < nx - 1; ++ii)
    array[ii * ny + jj] = 0;
}


void mpi_evolve(int pid, int np, double * array, int nx, int ny) {
  for (int ii = 1; ii < nx-1; ++ii) {
    for (int jj = 0; jj < ny; ++jj) {
      // check if boundary
      if (0 == pid && ii == 1)
        continue;
      if (pid == np-1 && ii == nx - 2)
        continue;
      if (jj == 0)
        continue;
      if (jj == ny - 1)
        continue;
      // evolve non boundary
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
  for (int ii = 1; ii < nx-1; ++ii) { // only real rows
    for (int jj = 0; jj < ny; ++jj) {
      std::cout << (ii-1 + pid*(nx-2)) * DELTA << " " << jj * DELTA << " " << array[ii * ny + jj] << "\n";
    }
    std::cout << "\n"; // separa las filas
  }
}



///////////////////// serial

void initial_conditions(Matrix &m) {
  for (int ii = 0; ii < N; ++ii) {
    for (int jj = 0; jj < N; ++jj) {
      m[ii * N + jj] = 1.0;
    }
  }
}

void boundary_conditions(Matrix &m) {
  int ii = 0, jj = 0;

  ii = 0;
  for (jj = 0; jj < N; ++jj)
    m[ii * N + jj] = 100;

  ii = N - 1;
  for (jj = 0; jj < N; ++jj)
    m[ii * N + jj] = 0;

  jj = 0;
  for (ii = 1; ii < N - 1; ++ii)
    m[ii * N + jj] = 0;

  jj = N - 1;
  for (ii = 1; ii < N - 1; ++ii)
    m[ii * N + jj] = 0;
}

void evolve(Matrix &m) {
  for (int ii = 0; ii < N; ++ii) {
    for (int jj = 0; jj < N; ++jj) {
      // check if boundary
      if (ii == 0)
        continue;
      if (ii == N - 1)
        continue;
      if (jj == 0)
        continue;
      if (jj == N - 1)
        continue;
      // evolve non boundary
      m[ii * N + jj] = (m[(ii + 1) * N + jj] + m[(ii - 1) * N + jj] +
                        m[ii * N + jj + 1] + m[ii * N + jj - 1]) /
                       4.0;
    }
  }
}

void print_gnuplot(const Matrix &m) {
  for (int ii = 0; ii < N; ++ii) {
    for (int jj = 0; jj < N; ++jj) {
      std::cout << ii * DELTA << " " << jj * DELTA << " " << m[ii * N + jj]
                << "\n";
    }
    std::cout << "\n";
  }
}


void init_gnuplot(void) {
  std::cout << "set contour " << std::endl;
  std::cout << "set terminal gif animate " << std::endl;
  std::cout << "set out 'anim.gif' " << std::endl;
}

void plot_gnuplot(const Matrix &m) {
  std::cout << "splot '-' w pm3d " << std::endl;
  print_gnuplot(m);
  std::cout << "e" << std::endl;
}

void print_matrix(const Matrix &m) {
  for (int ii = 0; ii < N; ++ii) {
    for (int jj = 0; jj < N; ++jj) {
      std::cout << m[ii * N + jj] << "  " ;
    }
    std::cout << "\n";
  }
}