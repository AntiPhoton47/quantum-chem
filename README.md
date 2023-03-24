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
If you find this program useful in your research, and you would like to acknowledge it in a journal publication, it is recommended to cite one of the following papers:
> P. A. LeMaitre and R. B. Thompson, ``Gaussian Basis Functions for an Orbital-Free-Related Density Functional Theory of Atoms''. Int. J. Quantum Chem. e27111 (2023).

> P. A. LeMaitre and R. B. Thompson, ``On the Origins of Spontaneous Spherical Symmetry-Breaking in Open-Shell Atoms Through Polymer Self-Consistent Field Theory''. J. Chem. Phys. 158 (6), 064301 (2023).

Here are example BibTeX entries:
```
@article{LeMaitre_Thompson2023,
  author  = {P. A. LeMaitre and R. B. Thompson},
  title   = {Gaussian Basis Functions for an Orbital-Free-Related Density Functional Theory of Atoms},
  journal = {Int. J. Quantum Chem.},
  pages = {e27111},
  year    = {(2023)},
  url  = {https://onlinelibrary.wiley.com/doi/full/10.1002/qua.27111},
}
@article{LeMaitre_Thompson2023ang,
  author  = {P. A. LeMaitre and R. B. Thompson},
  title   = {On the Origins of Spontaneous Spherical Symmetry-Breaking in Open-Shell Atoms Through Polymer Self-Consistent Field Theory},
  journal = {J. Chem. Phys.},
  volume = {158},
  number = {6},
  pages = {064301},
  year    = {(2023)},
  url  = {https://aip.scitation.org/doi/10.1063/5.0131364}
}
```
## Contributing
I am by no means a programming expert, so I welcome any suggestions for improving the program. 
