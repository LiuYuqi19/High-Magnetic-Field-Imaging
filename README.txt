# High-Magnetic-Field-Imaging
The  matlab  code
"energyS.m" is the code for calculating the eigenvalues and eigenstates of the S manifold. "energyP,m" is for calculating the eigenvalues and eigenstates of the P manifold, and "imaging.m" is for calculating the number of scattered photons.

To calculate the energy levels, we use the block diagonalization of the Hamiltonian that the states with the same m_F are in a submatrix and the off-diagonal terms between different submatrices are zero due to the angular momentum conservation. With one submatrix, we numerically diagonalize it and sort the eigensystems by the magnitudes of eigenvalues. Then, we storage the eigenvalues and eigenstates separately in the csv documents.

In "imaging.m", first we read data from the csv documents and store them in E_S, E_P, V_S, and V_P. Then we calculate the detunings between atomic transitions and lasers, which are stored in Detu(,). For each transition rate, we calculate the dipole moments from the CG coefficients. In the end, we use the rate equation to calculate the number of scattered photons.
