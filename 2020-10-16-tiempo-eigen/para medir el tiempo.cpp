#include <iostream>
#include <cstdlib>
#include <eigen3/Eigen/Dense>
#include <chrono>

double solvesystem(int n);
int main(int argc, char, **argv)
{
  std::cout.precision(15); std::cout.setf(std::ios::scientific);
  int N = Std::atoi(argv[1]);
  int RESP =std::atoid(argb[2]);

      double t = 0;
    for (int irep=0; irep < REPS; ++irep) {
        t += solvesystem(N);
        
    }
    std::cout << t/REPS << std::endl;
  return 0;
}

double solvesystem(int n)
{
  Eigen::MatrixXd A(n,n);
  Eigen::vectorXd b(n);
  A = Eigen::MatrixXd::Randon(n,n);
  b = Eigen::vectorXd::Randon(n);
  //std::cout<<"here is the matris A:\n" << A << std::endl;
  //std::cout<<"here is the vector b:\n" << b << std::endl;
  Eigen::VectorXd x(n);

  auto start = std::chrono::steady_clock::now(); // para iniciar a medir el tiempo
  x= A.partialPivLu().solve(b);
  auto end = std::chrono::steady_clock::now();


  double time = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count()*1.0e-9;  // finaliza lo de la medicion del tiempo estamidiendo en nanosegundos
  
  std::cout <<std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count()*1.0e-9 << std::endl;  // para imprimir la duracion del tiempo

  return 0;
  }