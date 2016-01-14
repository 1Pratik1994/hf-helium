# hf-helium
A C++ implementation of the Hartree-Fock solution of the ground state of the helium atom using Gaussian Type Orbitals (GTOs) as basis functions. The implementation uses LAPACK's LAPACKE_dsygv routine to solve the generalized eigenvalue problem. 

`hf-helium.cpp` can be compiled by invoking 
```
g++ -o hf-helium hf-helium.cpp -llapacke
```

