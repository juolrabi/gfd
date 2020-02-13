# Generalized finite differences (GFD) source code for solving large scale wave equations

This is a repository for a research project that aim to create an efficient numerical solver for large scale wave equations. The solver is based on the methodology called discrete exterior calculus (DEC). We call the outcoming method as generalized finite differences, because the solution method is closely related to finite-difference approaches.

The source code implements a library for the neccessary concepts and algorithms. User interface is still a mess. Sorry for that.

The development of the library has been active and in use since 2012 and is still under construction. There are several existing simulation results generated with this code or with previous versions of it. For more details, see [the projects web page](https://sites.google.com/jyu.fi/gfd).

The code in this repository is free software and licensed under the GNU General Public License version 3.

The programs are developed and run with Linux (Ubuntu). To build and run the code you need to install following packages:
- sudo apt install make
- sudo apt install g++
- sudo apt install libopenmpi-dev

The program samples can be built and run with following commands: 
- make: to build
- make run: to build and run
- make mpirun N: to run with N cores using MPI parallelization (not working with all samples)
- make clean: to clean build

Author of this project is Jukka Räbinä, who currently works in the 
- Faculty of Information Technology at the University of Jyväskylä.

The project is partially funded by
- Academy of Finland Grants No. 259925, 260076, and 295897,
- ERC Advanced Grant No. 320773.
