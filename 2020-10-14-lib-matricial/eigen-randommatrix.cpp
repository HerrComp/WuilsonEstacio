#include <iostream>
#include <eigen3/Eigen/Dense>

int main(int argc, char **argv)
{
    srand(1); // es para controlar la semilla 
    // unif entre -1 y 1
    Eigen::MatrixXd m = Eigen::MatrixXd::Random(3, 3);
    std::cout << m << std::endl;
    /* para pasar a una distribucion entre 0 y 1 tenemos que hacer los siguiente:
    Eigen::MatrixXd m = Eigen::MatrixXd::Random(3, 3);
    std::cout << m << std::endl;
     m=(m+1)/2
     std::cout << m << std::endl; 
     */

    /* para pasar entre a y b  
    double a=2.0, b= 3.0;
    m = a + (b-a)*m
    std::cout << m << std::endl;
     */ 

// ahora vamos a crear otra matriz de 2,2  osea cambiamos X por 2
    Eigen::Matrix2d m2;
    m2 << 1, 2,
        4, 4;
    std::cout << m2 << std::endl;  // para imprimir
    m2.setRandom();  // para convertir la matriz m2 en aleatoria distribuida entre -1 y 1
    std::cout << m2 << std::endl;

    return 0;
}