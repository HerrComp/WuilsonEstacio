#include <iostream>
#include <omp.h>

int main(int argc, char **argv)
{
  double x = 25; // global
#pragma omp parallel 
  {
    int tid = omp_get_thread_num();
    int nth = omp_get_num_threads();
    if (0 == tid) {
      std::cout << "nth: " << nth << std::endl;
    }
    double y = 56.8; // privada
    //std::cout << "Hello world, from tid: " << tid << std::endl;
    x = tid/34.5;
    //std::cout << "My tid:" << tid << std::endl;
    //std::cout << "My nth:" << nth << std::endl;
  }

  std::cout << x << std::endl;
  return 0;
}