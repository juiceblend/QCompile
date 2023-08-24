# QCompile
Compiles an arbitrary qubit unitary transformation into braids

Implementation of methods outlined in Kliuchnikov et al., [*Asymptotically Optimal Topological Quantum Compiling*](https://arxiv.org/pdf/1310.4150.pdf)

The main method is `compile_unitary` in `compile_unitary.py`, which takes in as input an arbitrary unitary matrix U and outputs a braid approximating that matrix to precision epsilon. 

An updated version (not yet fully functional for compile_unitary) is on branch `mpmath`.
