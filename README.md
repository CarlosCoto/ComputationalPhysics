# Computational Physics

Computational physics scripts repository

Numerical methods:

- Interpolation.ipynb: Lagrange interpolation with Fortran

- Integration.ipynb: Simpson's composite method with Fortran

- SystemOfEquations.ipynb: Resolution of system of equations by Gauss-Seidel Method with Fortran


The error at the bottom of the code blocks ("LFortran Exception: visit_Read() not implemented") is due to the lack of implementation of file IO in the runtime library of LFortran for Jupyter
https://fortran-lang.discourse.group/t/lfortran-exception-visit-read-not-implemented-error/1803
https://gitlab.com/lfortran/lfortran/-/issues/563

You can run these pieces of code on https://www.onlinegdb.com/ online Fortran compiler, or download a fortran compiler for your local system; most common is the the GNU Fortran compiler.
