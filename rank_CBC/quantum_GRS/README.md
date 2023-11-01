# quantum_cryptanalysis/rank_CBC/quantum_GRS

This repository provides the reference implementation of the computational costs for rank code-based cryptography against the quantum GRS algorithm. 

For example, if you run classical_GRS.c, then you execute as follows in this directory:
```
gcc -c classical_GRS.c -lgmp && gcc -o classical_GRS ../Rank_CBC_params.o basic_operation.o init_costs.o classical_GRS.o -lgmp && ./classical_GRS
```

Also, if you modify basic_operation.c or init_costs.c, then you execute as follows in this directory:
```
gcc -c basic_operation.c
```
or
```
gcc -c init_costs.c
```
