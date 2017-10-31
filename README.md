# TD-BEM #
A Parallel C++ Time-Domain Boundary Element code for elastodynamics analysis of solids;

### External libraries ###
* MPI for parallel programing;
* Petsc for solving large-scale linear equation system;
* Results can be written in VTK format and visualized with paraview;

### Key Boundary Element Techniques ###
* Regularization by Extracting the hyper-singular part of the Green's function;
* Singluar transformation for weakly-singluar element;
* Time marching with weighted collocation and higher-order projection;
