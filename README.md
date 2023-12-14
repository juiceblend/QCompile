# QCompile
Compiles an arbitrary qubit unitary transformation into braids

Implementation of methods outlined in Kliuchnikov et al., [*Asymptotically Optimal Topological Quantum Compiling*](https://arxiv.org/pdf/1310.4150.pdf)

The main methods are in `compile_rotation_cleaned.ipynb` and `compile_weave.ipynb`, which takes in as input an angle $\phi$ and outputs a braid/weave approximating $R_z(\phi)$ matrix to precision epsilon. 

Functions to manipulate the look of braids (e.g. to convert the printed output here to LaTeX) are in `Circuits.py`. The convention used in printing is that 1, 2 correspond to $\sigma_1, \sigma_2$ in the Kliuchnikov paper, while 3, 4 correspond to $\sigma_1^{-1}, \sigma_2^{-1}$, respectively.
