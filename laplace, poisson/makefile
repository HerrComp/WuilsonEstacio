
laplace.x: laplace-mpiv8.cpp
  mpic++ -std=c++17 -O3 laplace-mpiv8.cpp -o laplace.x

run: laplace.x runv8.sh
  mpic++ -std=c++17 -O3 laplace-mpiv8.cpp
  bash runv8.sh

datos.txt: laplace.mpi
  ./a.out 32 > datos.txt

clean:
  rm -f *.o *.mpi *.txt *.pdf