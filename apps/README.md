# ExasimApps

First run

```
/path/to/Exasim/build/text2code pdeapp.txt
```

Then run 

```
/path/to/Exasim/build/cput2cEXASIM pdeapp.txt                     (if you run on one CPU core)

/path/to/Exasim/build/gput2cEXASIM pdeapp.txt                     (if you run on one GPU)

mpirun -np $N /path/to/Exasim/build/cpumpit2cEXASIM pdeapp.txt    (if you run on many CPU cores)

mpirun -np $N /path/to/Exasim/build/gpumpit2cEXASIM pdeapp.txt    (if you run on many GPUs) 
```
