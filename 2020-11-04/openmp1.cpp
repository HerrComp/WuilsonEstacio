// se compila con: g++ -std=c++17 -fconcepts openmp1.cpp
#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <numeric>
#include <chrono>
//#include <eigen3/Eigen/Dense>

void print_elapsed(auto start, auto end);

int main(int argc, char, **argv) 
{
  std::cout.precision(15); std::cout.setf(std::ios::scientific);
  const int N = std::atoi(argv[1]);
  double * data = new double [N];

  auto star = std::chrono::steady_clock::now();
//fiil the aray
  for(int ii=0; ii < N; ++ii)
  {
    data[ii]=2*std::sin(ii)+std::log(ii+1);
  }
  auto end = std::chrono::steady_clock::now();
  //print_elapsed(start, end);

  double time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count()*1.0e-9;  // cuadra unidades de dedicion del tiempo aqui estamidiendo en nanosegundos

  //print
  std::cout << data[N/2] << std::endl;

delete [] data;
return 0;
}

void print_elapsed(auto start, auto end)
{
   std::cout << "Elapsed time in ns: " <<std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count()*1.0e-9 << std::endl;
}