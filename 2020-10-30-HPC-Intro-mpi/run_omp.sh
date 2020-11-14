echo "Running for all threads..."
for NTHREADS in 1 2 3 4 8 12 16; do echo  -n  "NTHREADS " OMP_NUM_THREADS=$NTHREADS;  ./a.out; done $> fulltimes.txt  

echo "Selecting elapased time.."
awk '{if (NR%2==1){gsub("elapsed"," ",$4); gsub("0:", " ", $4); print $1,$4}}' fulltimes.txt > time.txt

echo "Done"
TP1=$(awk '{if (NR==) print $2}' times.txt)

cat << EOF > plot.gp
    set term pdf
    set out "speedup.pdf"
    plot 'times.txt' u 1:($TP1/\$2) w lp 'speedup'
    set termpdf
    set out "efficiency.pdf"
    plot 'times.txt' u:1($TP1/\$2/\$1) w lp t 'parallel efficiency'
    EOF
    gnuplot plot gp
