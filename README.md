# FACUPD
MATLAB test code for direct factorization update of elliptic partial differential equations. This project is joint work with Dr. Jianlin Xia and Dr. Maarten V. de Hoop.

The folder FEM/ consists of functions to generate finite element matrices. All but localfem2d.m are selected from the nodal-dg package (https://github.com/tcew/nodal-dg). See the book Nodal Disontinuous Galerkin Methods by Jan S. Hesthaven and Tim Warburton for reference.

The folder SOL/ contains the factorization and solution functions.

An example is provided in main.m.
