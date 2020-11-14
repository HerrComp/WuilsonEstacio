for N in 1 2 10 100 500 1000 5000; do
	 echo -n "$N   "
	 ./eigen-Axb.x  $N 10  
done > datos.txt