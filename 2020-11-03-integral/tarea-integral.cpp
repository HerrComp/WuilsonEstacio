// for NT in 1 2 4 6 8; do echo "NTHREADS=$BT" export OMP_NUM_THREADS=$NT; time ./a.out 2; done
// 2> fulltimes.txt  // esta parte es si la quire guardar en un archivo txt
// g++ -std=c++17 -fconcepts -fopenmp tarea-integral.cpp

#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <numeric>
#include <chrono>

const int N = 1000000;
const double XMIN = 0.0;
const double XMAX = 10.0;
const double DX = (XMAX-XMIN)/N;

double f(double x);

void print_elapsed(auto start, auto end );

int main(int argc, char *argv[]) 
{
  std::cout.precision(15); 
  std::cout.setf(std::ios::scientific);
  
  int nthreads = std::atoi(argv[1]);

  double suma = 0.0;

#pragma omp parallel // for num_threads(nthreads) reduction(+:suma)  
  int nth = omp_get_num_threads(); //ask how many nthreads hay
  int tid = omp_get_thread_num(); // ask who i'm
  int SL = N/nthreads; // ask local size

  for (int ii = 0; ii < N; ++ii) {
    double x = XMIN + ii*DX;
    suma += DX*f(x); 
  }

  std::cout << "Integral calculada: " << suma << std::endl;
  std::cout << "Integral esperada : " << (XMAX*XMAX*XMAX - XMIN*XMIN*XMIN)/3 << std::endl;

  return 0;
}

double f(double x)
{
  return x*x;
}
