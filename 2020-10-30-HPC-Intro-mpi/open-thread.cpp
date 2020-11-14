#include <omp.h>
#include <stdio.h>


    int main(int argc, char *argv[]) {//PARAMETROS DESDE LA TERMINAL
  
      int nthreads, tid; //variables compartidas tid es quien sou yo 
  
      /* Fork a team of threads with each thread having a private tid variable */
    #pragma omp parallel private(tid) //private(tid) tid se esta declarando privada
      {
    
        /* Obtain and print thread id */
        tid = omp_get_thread_num();
        printf("Hello World from thread = %d\n", tid);
    
        /* Only master thread does this */
        if (tid == 0) 
        {
          nthreads = omp_get_num_threads();
          printf("Number of threads = %d\n", nthreads);
        }
    
      }  /* All threads join master thread and terminate */
  
    }