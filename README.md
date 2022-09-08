# Quantum Density Functional Theory Using Polymer Self-Consistent Field Theory
A Python program that computes the electronic densities and binding energies for isolated neutral atoms from an orbital-free density functional scheme based on polymer self-consistent field theory. The program uses a basis function expansion method to convert the model equations into a spectral form, where the basis functions are a product of radial Gaussians and real spherical harmonics. Further details about the algorithm, basis functions, and anything else can be found in the following thesis:
https://uwspace.uwaterloo.ca/handle/10012/18577

## Installation
The program uses some functions written in Fortran, which are wrapped using f2py and then imported as external modules. See the install file for instructions on how to set everything up on your machine.


Aditionally, the program requires the following external Python modules to be installed on your machine:
- NumPy
- SymPy
- SciPy
- Matplotlib
- PyVista

## Citing
If you find this program useful in your research, and you would like to acknowledge it in a journal publication, it is recommended to cite the following paper:
> P. A. LeMaitre and R. B. Thompson, ``Gaussian Basis Functions for a Polymer Self-Consistent Field Theory of Atoms''. arXiv:2208.09078 [physics.atom-ph] (2022).


Here is an example BibTeX entry:
```
@misc{LeMaitre_Thompson2022,
  author  = {P. A. LeMaitre and R. B. Thompson},
  title   = {Gaussian Basis Functions for a Polymer Self-Consistent Field Theory of Atoms},
  publisher = {arXiv},
  year    = {(2022)},
  url  = {https://arxiv.org/abs/2208.09078},
  doi = {10.48550/ARXIV.2208.09078}
}
```
## Contributing
I am by no means a programming expert, so I welcome any suggestions for improving the program. 
