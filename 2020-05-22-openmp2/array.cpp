// for NT in 1 2 3 4; do echo "NTHREADS=$BT" export OMP_NUM_THREADS=$NT; time ./a.out; done
// 2> fulltimes.txt  // esta parte es si la quire guardar en un archivo txt

// for NT in 1 2 3 4; do echo "NTHREADS=$BT" export OMP_NUM_THREADS=$NT; time ./a.out; done
// 2> fulltimes.txt  // esta parte es si la quire guardar en un archivo txt

#include <omp.h>
#include <iostream>
#include <cstdio>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cstdlib>


const int N = 1000000;

int main(int argc, char *argv[]) 
{
  std::cout.precision(15); 
  std::cout.setf(std::ios::scientific);

  std::vector<double> data(N);
  std::fill(data.begin(), data.end(), 1.5); 
  int nthreads = std::atoi(argv[1]); // para colocar en consola el numero de threads eg ./a.out 10
  std::vector<double> sumas(nthreads, 0.0);
  std::fill(sumas.begin(), sumas.end(), 0.0);

//inicia elcalculo de la particion de la memoria
#pragma omp parallel num_threads(nthreads)
  {
    int tid = omp_get_thread_num();
    int size = N/nthreads; // tamaño threads o particiones
    int imin = tid*size;
    for (int ii = imin; ii < imin + size; ++ii) {
      sumas[tid] += data[ii]; 
    }
  } //finaliza elcalculo de la particion de la memoria.

// para sumar todas las sumas totales de cada threads
  double suma = 0.0;
  for(auto x : sumas) {
    suma += x;
  }

// hace los caulculo totales
  std::cout << "mean = suma " << suma << std::endl;
  std::cout << "mean calculada: " << suma/N << std::endl;
  std::cout << "mean esperada : " << 1.5 << std::endl;
 
  return 0;
}
