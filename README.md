Welcome to MarLicS
------------------

 Marlics is a software written in C++ to to solve either the Beris-Edwards equation of nematodynamics without flow, or to find minimum energy states of the Landau-de Gennes free energy with the (FIRE) engine. Marlics comes with the most common setups where these simulations are performed, such as confined slab cells and spherical liquid crystal droplets. Different simulations can be performed passing a file describing the system geometry. FInally, we provide the most common anchoring energies as boundary conditions, such as Rapini-Papoular and Founier-Galatola.

The program takes as input a descriptive file giving the simulations parameters and initial conditions generating a series of different snapshots. Some examples of these input files can be found in the examples folder.

In the aux folder you can find some useful scripts with help visualize and start liquid crystal simulations. In renderize scripts you will find some scripts that renderives a view of the director in a given plane. There is also another script which generates iso surfaces of the thermodynamic order parameter in the same folder. In the boundaries folder the user can find some simple program, also written in C++, which can generate geometry files for the sphere and cylinder geometry. Both folders have a readme file explaining how the scripts/program are used.

File structure:
--------------

src/ -> folder containing MarLiCs the source code;

src/boundary.(cpp/h) -> a superclass used to define the anchoring energies;
src/boundary_%%%.(cpp/h) -> the class defining the anchoring type %%%;
src/container.(cpp/h) -> Responsible for printing the output CSV file;
src/driver.(cpp/h) -> Responsible for reading the input_file, setuping the simulation variables. Contains the definition of struct Simulation_Parameters, which contains all the require variables.
src/energy.(cpp/h) -> a superclass used to define the functional derivative to be used along with the integrator;
src/energy_ldg.(cpp/h) -> a class with the Landau de Gennes functional derivative;
src/geometry.(cpp/h) -> a superclass used to define the system geometry;
src/geometry_%%%.(cpp/h) -> the class which define the %%% geometry;
src/initial_conditions.(cpp/h) -> here we define all pre-programed initial conditions; 
src/integrator.(cpp/h) -> a superclass used to define the any temporal integrator
src/integrator_%%%.(cpp/h) -> the class which define the %%% integrator;
src/fire_integrator.(cpp/h) -> a class containing the FIRE method;
src/force.(cpp/h) -> a superclass used to define the force to be used along with the FIRE method;
src/force_ldg.(cpp/h) -> the class with the force defined after the LdG energy;

auxiliar/ -> some auxiliar files to aid the user;

auxiliar/renderize/ -> contain routines to renderize the directors and the isosurface of the parameter order S;
auxiliar/boundaries/ -> contain 2 examples off c++ files that can be used to produce boundary files;
examples/ -> a series of MarLiCs usage examples;
notebook/ -> Some Mathematica nb files used to evaluate symbolically the functional derivatives.
	
makefile -> a makefile used to compile the code
Windows_makefile -> a makefile to compile a windows executable under a GNU environment

MarlicsVS/ ->folder containing the files required to compile  marlics in the virtual studio environment.
