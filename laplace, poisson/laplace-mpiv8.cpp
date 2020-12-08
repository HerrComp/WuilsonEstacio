// mpic++ -std=c++17 -fsanitize=address -fconcepts -g -O3 laplace-mpivf.cpp
//-Werror -Wall
// mpic++ -std=c++17 -O3 laplace-mpiv8.cpp
// mpirun -np 4 ./a.out 12
//Wuilson Estacio y Yul villaba
/* gnuplot
set pm3d; set contour
set term pdf; set out 'matrix.pdf'
splot 'datos1.txt' w pm3d
*/

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include "mpi.h"

// constants
//const int N = int(L / DELTA) + 1;
int N = 10;
const double L = 1.479;
double DELTA = L/N;  // resolucion
const int STEPS = 200;



// MPI versions modela arreglo unidim
void print_matrix_slice(double * array, int nx, int ny);
void mpi_print_matrix(int pid, int np, double * array, int nx, int ny); // para generar el  arreglo de la matriz
void mpi_interchange_data(int pid, int np, double * array, int nx, int ny); // nx = Nl + 2
void mpi_initial_conditions(int pid, int np, double * array, int nx, int ny); // Valor inicial a la matrix
void mpi_boundary_conditions(int pid, int np, double * array, int nx, int ny); // propagacion de fronteras hacia el sistem
void mpi_evolve(int pid, int np, double * array, int nx, int ny); // para hacer la ecolucion de la matrix
void mpi_print_gnuplot(int pid, int np, double * array, int nx, int ny);
void mpi_print_gnuplot_slice(int pid, int np, double * array, int nx, int ny); // hacer animation

int main(int argc, char **argv) {
  
  N = std::atoi(argv[1]); // se sobre escribe
  DELTA = L/N; //se recalcual

  // iniciamos mpi
  int pid, np;
  MPI_Init(&argc, &argv); //Para que todos los procesos puedan leer la misma linea de comandos
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  // local array incluyendo los ghost
  int Nl = N/np; // np es el numero de procesos 
  double * data = new double [(Nl+2)*N] {0.0};
  // llenar con el pid actual
  // fill recibe varios argmument
  // solo necesita donde arranca y donde termina
  // esto es std::fill(inicia, finaliza, conque lo llenamos); que significa la memoria corrida sierto dato.
  std::fill(data, data + (Nl+2)*N, pid); 
  // pid es quien soy yo, np cuantos somos
  // initial and boundary conditions
  mpi_initial_conditions(pid, np, data, Nl+2, N);
  // Interchange data or infromation
  mpi_interchange_data(pid, np, data, Nl+2, N);
   // for de boundari conditions
  mpi_boundary_conditions(pid, np, data, Nl+2, N);

  mpi_interchange_data(pid, np, data, Nl+2, N);
  // evolucion en un loop
  for (int istep = 0; istep < STEPS; ++istep) {  
    mpi_evolve(pid, np, data, Nl+2, N);
    mpi_interchange_data(pid, np, data, Nl+2, N);
  }
  // imprimir
  mpi_print_matrix(pid, np, data, Nl + 2, N); // para poder imprimir la matrix
  mpi_print_gnuplot(pid, np, data, Nl+2, N);
  // para graficar en gnuplot
  delete [] data; // para liberar memoria
  MPI_Finalize();
  return 0;
}
// esto imprime la Matri completa
// yo imprimo mi matri
    //pido las matrices de los demas y las imprimo
    // de lo contrario envio mi matriz
void mpi_print_matrix(int pid, int np, double * array, int nx, int ny){
  if (0 == pid) {
    print_matrix_slice(array, nx, ny);
    double * buffer = new double [nx*ny];//para pedir memoria // buffer es un puntero . nx,ny es mi arreglo
    for (int ipid = 1; ipid < np; ++ipid) {
      MPI_Recv(buffer, nx*ny, MPI_DOUBLE, ipid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//recibir en buffer
      print_matrix_slice(buffer, nx, ny);//imprime en buffer
    }
    delete [] buffer; //liberar memoria
  } else {
    MPI_Send(array, nx*ny, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
}

// esto para cuadrar el arreglo de la matrix
void print_matrix_slice(double * array, int nx, int ny) { //esto recibe un puntero o el nombre del arreglo es  un puntero
  for (int ii = 0 ; ii < nx; ++ii) {
    // esta parte para saber cual fila es el ghost
    if (0 == ii or ii == nx-1) { 
      std::cout << "G : " ; 
    } else {
      std::cout << "  : " ;
    } // hasta aqui 
    for (int jj = 0 ; jj < ny; ++jj) {
      //std::cout << array[ii*ny + jj] << "  ";
      std::printf("%6.1lf ", array[ii*ny + jj]);
    }  // esto indica que imprima % numero de 6 casillas . 11 decimales, f de flotantes
    std::cout << "\n";
  }
}
/* para intercabiar informacion y enviar los datos hacia atras y hacia adelante 
    MPI_Send(donde arranca, cuantos elementos envia, tipo de datos, aquien se lo estoy enviando, el tag , MPI_COMM_WORLD);
    MPI_Recv(donde recibe, cuantos recive, tipo de datos, de quien recibes, el tag 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    */
void mpi_interchange_data(int pid, int np, double * array, int nx, int ny) { // nx = Nl + 2
  if (0 != pid) { //esto para que el cero no envia al anterior
    //enviar hacia atras
    MPI_Send(array + ny, ny, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD);
    MPI_Recv(array, ny, MPI_DOUBLE, pid-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  if (np-1 != pid) { // esto para que no envie el ultimo 
    MPI_Recv(array + ny*(nx - 1), ny, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(array + ny*(nx - 2), ny, MPI_DOUBLE, pid+1, 1, MPI_COMM_WORLD);
  }
}
// declaramos initial condition que recibe un arreglo
// llena todo de unos, depende de la condicion of frontera MPI
void mpi_initial_conditions(int pid, int np, double * array, int nx, int ny) {
  for (int ii = 1; ii < nx-1; ++ii) {
    for (int jj = 0; jj < ny; ++jj) {
      array[ii * ny + jj] = 1.0;
    }
  }
}
// condiciones de frontera con MPI 
// asocia valores con las condiciones de frontera 
/* 
(x,y), (x,L)=100, (L,y)=0, (0,y)=100, (x,0)=0,
*/
void mpi_boundary_conditions(int pid, int np, double * array, int nx, int ny) {
  int ii = 0, jj = 0;
  // (x,L)
  if (0 == pid) { 
    ii = 1; // fila real y empieza en 1 porque el 0 esta el ghost
    for (jj = 0; jj < ny; ++jj)
      array[ii * ny + jj] = 100;
  }
  // (x,0)
  if (pid == np-1) { // porque este pid el que tiene la ultima parte de la matriz
    ii = nx - 2; // fila real
    for (jj = 0; jj < ny; ++jj)
      array[ii * ny + jj] = 100;
  }
  // (0,y)  
  jj = 0;
  for (ii = 1; ii < nx - 1; ++ii)
    array[ii * ny + jj] = 0;
  // (L,y)
  jj = ny - 1;
  for (ii = 1; ii < nx - 1; ++ii)
    array[ii * ny + jj] = 0;
}

// recore la matrix y aplica algoritmos
// y verifica que no estemos modif la condition of frontera con MPI 
void mpi_evolve(int pid, int np, double * array, int nx, int ny) {
  for (int ii = 1; ii < nx-1; ++ii) { //qui estamos moviendonos desde el 1 hasta el nx-1 de las filas reales
    for (int jj = 0; jj < ny; ++jj) { // y aqui por todas las columnas
      // check if boundary y si estamos en el primer pid, primera parte de la matrix  y si estamos en la primera fila esa no la cambiamos porque es concidion de frontera.
      if (0 == pid && ii == 1) // esto nos dice que si soy el pid 0 me salto esa casilla
        continue;
      if (pid == np-1 && ii == nx - 2) // en esta parte estamos mirando la ultima parte de la matrix 
        continue;
      if (jj == 0) // si estoy en la primera columna la saltamos y eso aplica para allworld
        continue;
      if (jj == ny - 1) // si estoy in the last columna la saltamos y eso aplica para allword
        continue;
      // evolve non boundary si no estamos en una condicion de frontera se aplica el metodo de relajacion.
      array[ii * ny + jj] = (array[(ii + 1) * ny + jj] + array[(ii - 1) * ny + jj] +
			     array[ii * ny + jj + 1] + array[ii * ny + jj - 1]) / 4.0;
    }
  }
}

void mpi_print_gnuplot(int pid, int np, double * array, int nx, int ny)
{
  if (0 == pid) {
    mpi_print_gnuplot_slice(pid, np, array, nx, ny); //
    double * buffer = new double [nx*ny]; //para pedir memoria // buffer es un puntero . nx,ny es mi arreglo
    for (int ipid = 1; ipid < np; ++ipid) {
      MPI_Recv(buffer, nx*ny, MPI_DOUBLE, ipid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //recivi en buffer
      mpi_print_gnuplot_slice(ipid, np, buffer, nx, ny); //imprime en buffer
    }
    delete [] buffer; //liberar memoria
  } else {
    MPI_Send(array, nx*ny, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
}

void mpi_print_gnuplot_slice(int pid, int np, double * array, int nx, int ny) {
  for (int ii = 1; ii < nx-1; ++ii) { // esto tiene encuenta solo las filas reales
    for (int jj = 0; jj < ny; ++jj) { // aqui todas las columnas
      std::cout << (ii-1 + pid*(nx-2)) * DELTA << " " << jj * DELTA << " " << array[ii * ny + jj] << "\n"; // esto va imprimiendo nuestros arreglos
    }
    std::cout << "\n"; // esto es para separa las filas
  }
}



