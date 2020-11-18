#include <iostream>
#include <cstdlib>
#include <cmath>
#include <chrono>

void print_elapsed(auto start, auto end );

int main(int argc, char **argv)
{
  const int N = std::atoi(argv[1]);
  double * data = new double [N]; 

  auto start = std::chrono::steady_clock::now();
  // fill the array
  for (int ii = 0; ii < N; ++ii) {
    data[ii] = 2*std::sin(ii) + std::log(ii + 1);
  }
  auto end = std::chrono::steady_clock::now();
  print_elapsed(start, end);
  
  // print
  std::cout << data[N/2] << std::endl;
  
  delete [] data;
  return 0;
}

void print_elapsed(auto start, auto end )
{
  std::cout << "Elapsed time in ms: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << "\n";
}
std::cout << "Elapsed time in ns: " 
<<
std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count()*1.0e-9 << std::endl;