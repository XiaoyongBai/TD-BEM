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

### How to compile it ###
* Step 1: Specify the directories of include and library for MPI and PETSC in the CMakeFile
* Step 2: mkdir build; cd build; 
* Step 3: cmake ..; make

### How to use it ###
* mpirun -np 2 ./TD-BEM inputfile.txt
* Example of inputfile can be found in folder named "Geometry"

### How to do Finite Element - Boundary Element coupling ###
* Output the boundary element matrices from the main()
* Use the FE-BE code to do FE-BE coupling
