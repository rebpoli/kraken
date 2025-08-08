### What is this repository for? ###

* Phase Field Fracture Propagation 
* Integrated Phase field Advanced Crack Simulator (IPACS)
* http://users.ices.utexas.edu/~shlee/gallery_phasefield.html

### How do I get set up? ###

See http://dealii.org/8.4.1/readme.html for installation. 
= CMake version 2.8.8 or later
= deal.II-8.0 or higher
= Trilinos 11.4.1 or higher which is compatible with deal.II
= p4est 0.3.4.2
= open mpi  

### How to run ###

Setup ‘CMakeLists.txt’ file. 
( enables to switch  dim=2/dim=3 )
Need a parameter file to define the problem.
Need directories to save the output files         : output/
Need an appropriate mesh file in the directory    : input_msh/


### Who do I talk to? ###

This code is licensed under the "GNU GPL version 2 or later” and deal.ii. See https://www.gnu.org/licenses/gpl-2.0.html and https://www.dealii.org/license.html. Copyright 2013-2014-2015-2016: Thomas Wick and Timo Heister and Sanghyun Lee. This code was modified by Sanghyun Lee 2014-2015. This code was modified by Sanghyun Lee and Thomas Wick 2016.

IPACS is a research code and general parameter specifications have not been included and require code modifications. 
Recently, it is used by Ph.D students Sogo S. and Mohammad J.