#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
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

void print_matrix_slice(double * array, int nx, int ny);
void mpi_print_matrix(int pid, int np, double * array, int nx, int ny);
void mpi_interchange_data(int pid, int np, double * array, int nx, int ny); // nx = Nl + 2

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
  // communication
  mpi_interchange_data(pid, np, data, Nl+2, N);
  // imprimir
  mpi_print_matrix(pid, np, data, Nl + 2, N);

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
      std::cout << array[ii*ny + jj] << "  ";
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
    MPI_Send(array + ny*(nx - 2), ny, MPI_DOUBLE, pid+1, 1, MPI_COMM_WORLD);
    MPI_Recv(array + ny*(nx - 1), ny, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}






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