// for NT in 1 2 3 4; do echo "NTHREADS=$BT" export OMP_NUM_THREADS=$NT; time ./a.out; done
// 2> fulltimes.txt  // esta parte es si la quire guardar en un archivo txt

#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <cmath>

const int N = 100000000;
const double XMIN = 0.0;
const double XMAX = 1000.0;
const double DX = (XMAX-XMIN)/N;

double f(double x);

int main(int argc, char *argv[]) 
{
  std::cout.precision(15); 
  std::cout.setf(std::ios::scientific);
  
  int nthreads = std::atoi(argv[1]);

  double suma = 0.0;
// <<<<<<< HEAD
=======
#pragma omp parallel for num_threads(nthreads) reduction(+:suma)  
//>>>>>>> 5be1a175fae6e757d3c6936396a3b8a49742197b
  for (int ii = 0; ii < N; ++ii) {
    double x = XMIN + ii*DX;
    suma += DX*f(x); 
  }
// <<<<<<< HEAD
  std::cout << "Integral calculada: " << suma << std::endl;
  std::cout << "Integral esperada : " << (XMAX*XMAX*XMAX - XMIN*XMIN*XMIN)/3 << std::endl;
=======
  
  std::cout << "Integral calculada: " << suma << std::endl;
  //std::cout << "Integral esperada : " << (XMAX*XMAX*XMAX - XMIN*XMIN*XMIN)/3 << std::endl;
//>>>>>>> 5be1a175fae6e757d3c6936396a3b8a49742197b

  return 0;
}

double f(double x)
{
<<<<<<< HEAD
  return x*x;
=======
  return x*std::sin(x)/(2+std::cos(x));
//>>>>>>> 5be1a175fae6e757d3c6936396a3b8a49742197b
}