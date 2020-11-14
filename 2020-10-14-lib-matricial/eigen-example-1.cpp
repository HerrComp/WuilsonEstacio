#include <iostream>
#include <eigen3/Eigen/Dense>

int main(int argc, char **argv)
{
  // el X indica que es dinamica  osea que la  matrix puede cambiar de tamaño y el d  que es de presicion double
    Eigen::MatrixXd m(2, 2); 
    m(0, 0) = 3;
    m(1, 0) = 2.5;
    m(0, 1) = -1;
    m(1, 1) = m(1, 0) + m(0, 1);
    std::cout << m << std::endl; // para imprimir 
    
    return 0;
}