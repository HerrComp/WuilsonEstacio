for NT in 1 4 6 16 22; do echo "NTHREADS=$BT" export OMP_NUM_THREADS=$NT; time ./a.out 1000000 10; done > NT.txt;
