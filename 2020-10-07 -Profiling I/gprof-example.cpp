// 

#include <cstdio>
#include <cstdlib>
#include <chrono> // para poder medir el tiempo
#include <iostream>

/*
   Tests cache misses.
*/

void print_elapsed(auto start, auto end );
void option1(long * data, int nx, int ny);
void option2(long * data, int nx, int ny);
void option3(long * data, int nx, int ny);

int main(int argc, char **argv)
{
  if (argc < 3){
    printf("Usage: cacheTest sizeI sizeJ\nIn first input.\n");
    return 1;
  }
  long sI = atoi(argv[1]);
  long sJ = atoi(argv[2]);

  
  printf("Operating on matrix of size %ld by %ld\n", sI, sJ);

  auto start = std::chrono::steady_clock::now(); // DESDE QUI EMPIESA  A MEDIR EL RELOG
  long *arr = new long[sI*sJ]; // double array.
  auto end = std::chrono::steady_clock::now(); // A QUI FINALIZA LA MEDICION DEL RELOG
  print_elapsed(start, end);
  
  // option 1
  start = std::chrono::steady_clock::now();
  option1(arr, sI, sJ); // aqui se coloca a lo que se le va a medir el timpo de ejecucion
  end = std::chrono::steady_clock::now();
  print_elapsed(start, end);
  
  // option 2
  start = std::chrono::steady_clock::now();
  option2(arr, sI, sJ);
  end = std::chrono::steady_clock::now();
  print_elapsed(start, end);

  // option 3
  start = std::chrono::steady_clock::now();
  option3(arr, sI, sJ);
  end = std::chrono::steady_clock::now();
  print_elapsed(start, end);

  printf("%ld\n", arr[0]);
  
  return 0;
}
  

void print_elapsed(auto start, auto end )
{
  std::cout << "Elapsed time in ms: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() // paraimprimir laduracion del tiempo en milisegundos
            << "\n";
}

void option1(long * data, int nx, int ny) // accede primero a columnas luego a filas.
{
  for (long i=0; i < nx; ++i)
    for (long j=0; j < ny; ++j)
      data[(i * (ny)) + j ] = i;
}

void option2(long * data, int nx, int ny) // accede primero a filas luego a colmnas
{
  for (long i=0; i < nx; ++i)
    for (long j=0; j < ny; ++j)
      data[(j * (nx)) + i ] = i;
}

void option3(long * data, int nx, int ny) // aqui se hace como un arreglo adimensional. como esta gasta menor tiempo es lamejor opcion para hacer el codigo.
{
  for (int i=0; i < nx*ny; ++i) data[i] = i; 
}