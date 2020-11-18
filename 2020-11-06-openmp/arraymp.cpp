//compila con: g++ -std=c++17 -fconcepts -fopenmp arraymp.cpp

// se corre con: time ./a.out 1000000 NTH; done

// for NT in 1 4 6 16 22; do echo "NTHREADS=$BT" export OMP_NUM_THREADS=$NT; time ./a.out 1000000 #; done

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <omp.h>

void print_elapsed(auto start, auto end );

int main(int argc, char **argv)
{
  std::cout.precision(15); std::cout.setf(std::ios::scientific);
  const int N = std::atoi(argv[1]);
  //const int NTH = std::atoi(argv[2]);
  double * data = new double [N]; 

//begin time 
  auto start = std::chrono::steady_clock::now();
  // fill the array
  #pragma omp parallel // num_threads(NTH)
  {
    int nth = omp_get_num_threads(); //ask how many nthreads hay
    int tid = omp_get_thread_num(); // ask who i'm
    int SL = N/nth; // ask local size
    
    for (int ii = tid*SL; ii < tid*SL + SL; ++ii) {
      data[ii] = 2*std::sin(ii) + std::log(ii + 1);
    }
  }
  auto end = std::chrono::steady_clock::now();
  print_elapsed(start, end); //end time
  
  // print
  std::cout << data[N/2] << std::endl;
  
  delete [] data;
  return 0;
}

void print_elapsed(auto start, auto end )
{
  //print time
  std::cout << "Elapsed time in ms: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << "\n";
}