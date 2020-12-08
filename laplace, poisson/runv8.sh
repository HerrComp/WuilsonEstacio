# bash runv8.sh
for N in 1 2 3 4 6 8 10 ; do
    echo -n "$N  "
    time ./a.out 40
done