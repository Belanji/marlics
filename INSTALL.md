
		COMPILING MARLICS - LINUX

In order to compile the MarLiCs program, the user will need the GSL (GNU Scientific Library) and the OpenMP library installed in a linux machine.
MarLiCs can be compiled using the comand 'make' or 'make marlics' in this folder in order to compile a executable named marlics, using g++.
If the target machine, the one execunting the simulation, processor differs from the native machine, the one compiling the code, we recomend to remove the "-march=native" argument.
We also provide a template, as comented lines, in the makefile in order to compile the MarLiCs using the icpc compiler.

MarLiCs was programed to be executed in linux based OS but, if necessary, this program can be executed inside of a Windows machine.

		COMPILING MARLICS - LINUX TO WINDOWS

MarLiCs can be compiled in a Linux Machine to be executed in a Windows one using the MingW compiler, which can be installed via apt-get in Ubuntu or Debian.
After installing MingW, be sure to compile the GSL to the MinGW as well. This can be done using the following steps:

Download the gsl library in https://www.gnu.org/software/gsl/ and extract.
In the terminal go to extracted folder and create and move to a build directory using "mkdir build; cd build".
Prepare the instalation using:
	./configure --prefix=/usr/x86_64-w64-mingw32 --host=mingw64
and compile the library using:
	make && sudo make install

After these steps, MingW should be ready to compile a Windows version of MarLiCs.
The user can use the makefile_Windows, with the comand "make -f makefile_Windows" 

		COMPILING MARLICS - WINDOWS TO WINDOWS

MaLiCs can also be compiled directly in a Windows, but, to do so the user need to install the GSL library in the Windows enviroment.
One way to do this is using the vcpkg, which can be done with the following steps:

With the Visual Studio installed, open the PowerShell go to a suitable directory, says D:, to install the vcpkg and the GSL then use the commands:

	mkdir dev
	cd dev

and clone the vcpkg repository with:

	git clone https://github.com/microsoft/vcpkg.git
	
with the repository in hands do:

	cd vcpkg
	.\bootstrap-vcpkg.bat
	.\vcpkg.exe integrate install

to setup the vcpkg and then:

	.\vcpkg.exe install gsl gsl:x64-windows

in order to install the gsl.
After that, start a new session in the Visual Studio and open the project MarlicsVS.sln in the MarlicsVS folder of the repository.
NOTE: Until July, 2021, Visual Studio can't produce a parallel version of the Marlics, since MarLiCs needs the openMP 4.5 or above.



