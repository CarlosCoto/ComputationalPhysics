# Computational Physics

Computational physics scripts repository

Numerical methods (with Fortran):

- Interpolation.ipynb: Lagrange interpolation

- Integration.ipynb: Simpson's composite method

- SystemOfEquations.ipynb: Resolution of system of equations by Gauss-Seidel method

- ZeroesOfFunctions.ipynb : Search for zeroes of functions using the Secant method

- Overflow.ipynb: determine overflow in simple and double precision

- ssolar.f: Program that simulates the Solar System using Verlet algorithm for solving differential equations


The error at the bottom of the code blocks ("LFortran Exception: visit_Read() not implemented") is due to the lack of implementation of file IO in the runtime library of LFortran for Jupyter
https://fortran-lang.discourse.group/t/lfortran-exception-visit-read-not-implemented-error/1803
https://gitlab.com/lfortran/lfortran/-/issues/563

You can run these pieces of code on https://www.onlinegdb.com/ online Fortran compiler, or download one for your local system; most common is the GNU Fortran compiler.
