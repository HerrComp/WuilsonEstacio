// for NT in 1 2 3 4; do echo "NTHREADS=$BT" export OMP_NUM_THREADS=$NT; time ./a.out; done
// 2> fulltimes.txt  // esta parte es si la quire guardar en un archivo txt

//awk '{if (NR%2==1){gsub("elapsed"," ",$3); print $3}}' fulltimes.txt > time.txt

//bash run_omp.sh // se debe crear un archivo omp.sh  con todo lo anterior para ejecutar bash

#include <stdio.h>
#include <omp.h>
#include <iostream>
#include <cmath>

int main(void)
{
  const int N =8000000;
  int * data = nullptr;
  data = new int [N];

#pragma omp parallel for //para que corra en paralelo
  for(int i = 0; i < N; ++i){
    data[i] = 2*std::sin(i)*sin(i);
  }

  std::cout << data[N/2] <<std::endl;

  delete [] data;

  return 0;
}