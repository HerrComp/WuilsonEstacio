//mpic++ -std=c++17 -fsanitize=address -fconcepts -g -o3 laplace-v1.cpp
//-Werror -Wall

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include "mpi.h"

// constants
//const int N = int(L / DELTA) + 1;
const int N = 10;
const double L = 1.479;
const double DELTA = L/N; // resolucion
const int STEPS = 200;

typedef std::vector<double> Matrix;

void initial_conditions(Matrix &m); 
void boundary_conditions(Matrix &m);
void evolve(Matrix &m);// propagacion
void print_gnuplot(const Matrix &m);
void print_matrix(const Matrix &m);
void init_gnuplot(void); //impr commandos
void plot_gnuplot(const Matrix &m);

//
int main(void) {
  Matrix data(N * N);//se crea el arreglo
  initial_conditions(data);
  boundary_conditions(data);

  //init_gnuplot(); // para simulacion
  for (int istep = 0; istep < STEPS; ++istep) {
    evolve(data);
    //plot_gnuplot(data);
  }
  //print_gnuplot(data);
  print_matrix(data);
  
  return 0;
}
// llena todo de unos
void initial_conditions(Matrix &m) {
  for (int ii = 0; ii < N; ++ii) {
    for (int jj = 0; jj < N; ++jj) {
      m[ii * N + jj] = 1.0;
    }
  }
}
// asocia valores en los bordes 
/* 
(x,y), (0,y)=0,(l,y)=0,(x,0)=0,(x,L)=100
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

}// recore la matrix y app algorit
void evolve(Matrix &m) {
  for (int ii = 0; ii < N; ++ii) {
    for (int jj = 0; jj < N; ++jj) {
      // check if boundary
      if (ii == 0)
        continue; //
      if (ii == N - 1)
        continue;
      if (jj == 0)
        continue;
      if (jj == N - 1)
        continue;
      // evolve non boundary,si no estamos en una condicion de frontera se aplica elmetodo de relajacion
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