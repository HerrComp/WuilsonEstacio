#include <iostream>
#include <cstdlib>
#include <eigen3/Eigen/Dense>

int main(int argc, char **argv)
{
   
    Eigen::Matrix3d A; // esto es una matris de 3,3  porque X=3
    Eigen::Vector3d b;
    A = << 1,2,3  4,5,6,  7,8,10;
    b = << 3, 3, 4;
    //std::cout << "Here is the matrix A:\n" << A << std::endl;
    //std::cout << "Here is the vector b:\n" << b << std::endl;
    Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);
    //std::cout << "The solution is:\n" << x << std::endl;
    std::cout << "Error: " << (A*x - b).norm() << std::endl;

}