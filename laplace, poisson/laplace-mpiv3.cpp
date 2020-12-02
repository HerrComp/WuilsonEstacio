// mpic++ -std=c++17 -fsanitize=address -fconcepts -g -O3 laplace-mpiv3.cpp
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
const double L = 1.500;
double DELTA = L/N; // resolucion
const int STEPS = 200;

typedef std::vector<double> Matrix; //modela arreglo unidim

void initial_conditions(Matrix &m); // Valor inicial a la matrisx
void boundary_conditions(Matrix &m);
void evolve(Matrix &m); // propagacion de fronteras hacia el sistem
void print_gnuplot(const Matrix &m);
void print_matrix(const Matrix &m);
void init_gnuplot(void); //impr commandos
void plot_gnuplot(const Matrix &m); // hacer enimation

void print_matrix_slice(double * array, int nx, int ny);
void mpi_print_matrix(int pid, int np, double * array, int nx, int ny);
void mpi_interchange_data(int pid, int np, double * array, int nx, int ny); // nx = Nl + 2

int main(int argc, char **argv) {
  
  N = std::atoi(argv[1]);
  DELTA = L/N;

  int pid, np;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  // local array
  int Nl = N/np;
  double * data = new double [(Nl+2)*N] {0.0};
  // llenar con el pid actual
  // fill recibe varios argmument
  // solo necesita donde arranca y donde termina
  // std::fill(inicia, finaliza, conque lo llenamos); 
  std::fill(data, data + (Nl+2)*N, pid); 
  // communication
  mpi_interchange_data(pid, np, data, Nl+2, N);
  // imprimir //pid es quien soy yo, 
  // np cuantos somos
  mpi_print_matrix(pid, np, data, Nl + 2, N);

  delete [] data; // para liberar memoria
  MPI_Finalize();
  return 0;
}
// esto imprime la Matri completa
// yo imprimo mi matri
    //pido las matrices de los demas y las imprimo
    // de lo contrario envio mi matriz
void mpi_print_matrix(int pid, int np, double * array, int nx, int ny){
  if (0 == pid) {
    print_matrix_slice(array, nx, ny);
    double * buffer = new double [nx*ny]; //para pedir memoria // buffer es un puntero . nx,ny es mi arreglo
    for (int ipid = 1; ipid < np; ++ipid) {
      MPI_Recv(buffer, nx*ny, MPI_DOUBLE, ipid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //recivi en buffer
      print_matrix_slice(buffer, nx, ny); //impr en buffer
    }
    delete [] buffer; //liberar memoria
  } else {
  //} else {
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

// llena todo de unos, depende de la condicion of frontera
void initial_conditions(Matrix &m) {
  for (int ii = 0; ii < N; ++ii) {
    for (int jj = 0; jj < N; ++jj) {
      m[ii * N + jj] = 1.0;
    }
  }
}
// asocia valores en los bordes 
/* 
(x,y), (0,y)=0,(L,y)=0,(x,0)=0,(x,L)=100
*/
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
// recore la matrix y app algorit
//y verifica que no estemosmodif la condition of frontera
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
      // evolve non boundary, si no estamos en una condicion de frontera se aplica elmetodo de relajacion
      m[ii * N + jj] = (m[(ii + 1) * N + jj] + m[(ii - 1) * N + jj] +
                        m[ii * N + jj + 1] + m[ii * N + jj - 1]) /
                       4.0; //que es que el valor de mi selda es
    } //el promedio delos valores de los vecionos 
  }
}

void print_gnuplot(const Matrix &m) {
  for (int ii = 0; ii < N; ++ii) {
    for (int jj = 0; jj < N; ++jj) {
      std::cout << ii * DELTA << " " << jj * DELTA << " " << m[ii * N + jj] //imprime condenada en i,j y en valor de la Mtx hay
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
