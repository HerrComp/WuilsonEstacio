#include <stdio.h>
#include <omp.h>
#include <iostream>

int main(void)
{
  const int N =1000000;
  int data[N]={0};

#pragma omp parallel for //para que corra en paralelo
  for(int i = 0; i < N; ++i){
    data[i]=2*i;
  }

  std::cout << data[N/2] <<std::endl;

  return 0;
}