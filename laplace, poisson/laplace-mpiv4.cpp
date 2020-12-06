// mpic++ -std=c++17 -fsanitize=address -fconcepts -g -O3 laplace-mpiv4.cpp
//-Werror -Wall
// mpirun -np 4 ./a.out 10

//mpic++ -std=c++17 -O3 laplace-mpiv4.cpp

//Wuilson Estacio
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "mpi.h"
#include <cstdio>
//#include <chrono>


// constants
//const int N = int(L / DELTA) + 1;
int N = 10;
const double L = 1.500;
double DELTA = L/N; // resolucion
const int STEPS = 200;

typedef std::vector<double> Matrix;
//modela arreglo unidim

void initial_conditions(Matrix &m);// Valor inicial a la matrix
void boundary_conditions(Matrix &m);
void evolve(Matrix &m); // propagacion de fronteras hacia el sistem
void print_gnuplot(const Matrix &m);
void print_matrix(const Matrix &m);
void init_gnuplot(void); //impr commandos
void plot_gnuplot(const Matrix &m); // hacer enimation


// For versions MPI
void print_matrix_slice(double * array, int nx, int ny);
void mpi_print_matrix(int pid, int np, double * array, int nx, int ny);
void mpi_interchange_data(int pid, int np, double * array, int nx, int ny); // nx = Nl + 
void mpi_initial_conditions(int pid, int np, double * array, int nx, int ny);
void mpi_boundary_conditions(int pid, int np, double * array, bool * dat,  int nx, int ny);
void mpi_evolve(int pid, int np, double * array, bool * dat, int nx, int ny);


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
  // initial conditions
  mpi_initial_conditions(pid, np, data, Nl+2, N);
  // Interchange data or infromation
  mpi_interchange_data(pid, np, data, Nl+2, N);
  // imprimir //pid es quien soy yo, 
  // np cuantos somos
  mpi_print_matrix(pid, np, data, Nl + 2, N);

  delete [] data;
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
    double * buffer = new double [nx*ny];//para pedir memoria // buffer es un puntero . nx,ny es mi arreglo
    for (int ipid = 1; ipid < np; ++ipid) {
      MPI_Recv(buffer, nx*ny, MPI_DOUBLE, ipid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//recivi en buffer
      print_matrix_slice(buffer, nx, ny);//impr en buffer
    }
    delete [] buffer;  //liberar memoria
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

/* para intercabiar informacion y enviar los datos hacia atras y hacia adelante 
    MPI_Send(donde arranca, cuantos elementos envia, tipo de datos, aquien se lo estoy enviando, el tag , MPI_COMM_WORLD);
    MPI_Recv(donde recibe, cuantos recive, tipo de datos, de quien recibes, el tag 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    */
void mpi_interchange_data(int pid, int np, double * array, int nx, int ny) { // nx = Nl + 2
  if (0 != pid) { //esto para que el cero no envia al anterior
    MPI_Send(array + ny, ny, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD);
    MPI_Recv(array, ny, MPI_DOUBLE, pid-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  if (np-1 != pid) {// esto para que no envie el ultimo 
    MPI_Send(array + ny*(nx - 2), ny, MPI_DOUBLE, pid+1, 1, MPI_COMM_WORLD);
    MPI_Recv(array + ny*(nx - 1), ny, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}

// declaramos initial condition que recibe un arreglo
// llena todo de unos, depende de la condicion of frontera MPI
void mpi_initial_conditions(int pid, int np, double * array, int nx, int ny) {
  for (int ii = 1; ii < nx-1; ++ii) {
    for (int jj = 0; jj < ny; ++jj) {
      array[ii * ny + jj] = 1.0;
    }
  }
}

// declaramos initial condition que recibe un arreglo
// llena todo de unos, depende de la condicion of frontera 
void initial_conditions(Matrix &m) {
  for (int ii = 0; ii < N; ++ii) {
    for (int jj = 0; jj < N; ++jj) {
      m[ii * N + jj] = 1.0;
    }
  }
}

//condiciones de frontera con MPI 
// asocia valores con las condiciones fe frontera 
/* 
(x,y), (0,y)=0,(L,y)=0,(x,0)=0,(x,L)=100
*/
void mpi_boundary_conditions(int pid, int np, double * array, bool * dat,  int nx, int ny) {
  int ii = 0, jj = 0;

  if (0 == pid) {
    ii = 1; 
    for (jj = 0; jj < ny; ++jj)
      array[ii * ny + jj] = 100;
      dat[ii*ny + jj] = true;
  }

  if (pid == np-1) {
    ii = nx - 2; 
    for (jj = 0; jj < ny; ++jj)
      array[ii * ny + jj] = 0;
      dat[ii*ny + jj] = true;
  }
  
  jj = 0;
  for (ii = 1; ii < nx - 1; ++ii){
    array[ii * ny + jj] = 0;
    dat[ii*ny + jj] = true;
  }

  jj = ny - 1;
  for (ii = 1; ii < nx - 1; ++ii)
    array[ii * ny + jj] = 0;
    dat[ii*ny + jj] = true;
}

/*
void boundary_conditions(Matrix &m) {
  int ii = 0, jj = 0;

  ii = 0;
  for (jj = 0; jj < N; ++jj)
    m[ii * N + jj] = 100.0;

  ii = N - 1;
  for (jj = 0; jj < N; ++jj)
    m[ii * N + jj] = 0.0;

  jj = 0;
  for (ii = 1; ii < N - 1; ++ii)
    m[ii * N + jj] = 100.0;

  jj = N - 1;
  for (ii = 1; ii < N - 1; ++ii)
    m[ii * N + jj] = 0.0;
}
*/

// recore la matrix y app algorit
//y verifica que no estemos modif la condition of frontera con MPI 
void mpi_evolve(int pid, int np, double * array, bool * dat, int nx, int ny) {
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
      if(dat[ii * ny + jj] == true)
        continue;
      // evolve non boundary si no estamos en una condicion de frontera se aplica elmetodo de relajacion.
      array[ii * ny + jj] = (array[(ii + 1) * ny + jj] + array[(ii - 1) * ny + jj] +
			     array[ii * ny + jj + 1] + array[ii * ny + jj - 1]) / 4.0;
    }
  }
}
/*
//recore la matrix y app algorit
//y verifica que no estemos modif la condition of frontera
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
                       4.0;
    }
  }
}
*/

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