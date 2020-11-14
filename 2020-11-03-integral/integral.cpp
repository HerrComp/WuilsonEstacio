// for NT in 1 2 3 4; do echo "NTHREADS=$BT" export OMP_NUM_THREADS=$NT; time ./a.out; done
// 2> fulltimes.txt  // esta parte es si la quire guardar en un archivo txt

#include "omp.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <numeric>


const int N = 100000000;
const double XMIN = 0.0;
const double XMAX = 1000.0;
const double DX = (XMAX-XMIN)/N;

double f(double x);

int main(void) 
{
  std::cout.precision(15); 
  std::cout.setf(std::ios::scientific);
  
  double suma = 0.0;
  std::vector<double> partial;

//#pragma omp parallel for num_threads(nthreads) reduction(+:suma)  
#pragma omp parallel
{
  int tid = omp_get_thread_num();
  int nthreads = omp_get_thread_num();

  if(0==tid){
    partial.resize(nthreads);
  }
  #pragma omp barrier

  int size = N/nthreads;
  int init = tid*size;

  for (int ii = init; ii < init + size; ++ii) {
    double x = XMIN + ii*DX;
    partial[tid] += f(x);
  }
  partial[tid] *= DX;
}

std::cout << std::accumulate(partial.begin(), partial.end(), 0.0) << "\n";

//std::cout << "Integral esperada : " << (XMAX*XMAX*XMAX - XMIN*XMIN*XMIN)/3 << std::endl;

return 0;
}

double f(double x)
{
  return x*x;

  //return x*std::sin(x)/(2+std::cos(x));

}
