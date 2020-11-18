// for NT in 1 2 3 4; do echo "NTHREADS=$BT" export OMP_NUM_THREADS=$NT; time ./a.out; done
// 2> fulltimes.txt  // esta parte es si la quire guardar en un archivo txt

// for NT in 1 2 3 4; do echo "NTHREADS=$BT" export OMP_NUM_THREADS=$NT; time ./a.out; done
// 2> fulltimes.txt  // esta parte es si la quire guardar en un archivo txt

#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <numeric>


int main(int argc, char *argv[]) 
{
  std::cout.precision(15); 
  std::cout.setf(std::ios::scientific);

  int nthreads, tid;

  std::cout << "(0) Addres of tid:" << &tid << std::endl; // imprime la direcion de variable de una variable


#pragma omp parallel private(tid)
{
  //for num_threads(nthreads) reduction(+:suma)  
  tid = omp_get_thread_num();
  std::printf("Hello word from thread = %d\n",tid);

  if (tid == 0)
  {
    nthreads = omp_get_num_threads();
    std::printf("Numbe of threads =%d\n", nthreads);
  }
  
  std::cout << "Addres of tid:" << &tid << std::endl; // imprime la direcion de variable de una variable

}

nthreads = omp_get_num_threads();
std::printf("Numbe of threads =%d\n", nthreads);

}

