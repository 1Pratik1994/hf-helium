# hf-helium
A C++ implementation of the Hartree-Fock solution of the ground state of the helium atom using Gaussian Type Orbitals (GTOs) as basis functions <sup>1</sup>. The implementation uses LAPACK's LAPACKE_dsygv routine to solve the generalized eigenvalue problem

```
F * d = E * S * d
```
where `F` is the Fock matrix, `S` the overlap matrix, `d` the coefficients vector and `E` the energy eigenvalues. The experimental value for the ground state energy of the helium atom is `E0 = −2.903 hartree`; `hf-helium` returns the approximated value `E0 = −2.855160356 hartree`.

`hf-helium.cpp` can be compiled by invoking 
```
g++ -o hf-helium hf-helium.cpp -llapacke
```
<sup>1</sup>: The implementation is based on the information presented as part of Exercise 9 of Professor Matthias Troyer's *Computational Quantum Physics* lecture (spring semester 2015, ETH Zürich)
