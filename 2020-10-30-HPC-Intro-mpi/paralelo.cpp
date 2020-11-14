   // OMP_NUM_THREADS=4 ./a.out
#include <stdio.h>
#include <omp.h>
    
    int main(void)
    {
      int a = 9;
    #pragma omp parallel private (a) //num_threads (4)
    {
      double x = 0;
      printf("Hello, world.\n");
      a += 3;
      printf("%p\n", &x); //para generar impresionces de memoria
      printf("%p\n", a);

    }
      printf("%p\n", a);
      return 0;
    }