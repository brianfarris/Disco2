Disco2
======

Build Instructions

1) Copy `Makefile.in.template` into a new file `Makefile.in`.  Modify `Makefile.in` to accomodate your machine's layout, currently this only requires pointing it to your (parallel) HDF5 installation and setting the flags for performance or development.

2) Type `make` to compile the binary.

---

Run Instructions

The disco binary lives in `bin/`.  It takes a single command-line option: a parameter file.  Example parameter files are found in `parfiles/`.  If you're in the top directory you can run the vortex example simply as:
```
$ bin/disco parfiles/code_tests/vortex.par
```
Output files will be written into the directory where you executed disco.

---

Notes

* The overall approach here is to mimic the sort of object-oriented design philosophy that you would typically get in C++, using only C. The way we do this is by organizing data and functions into categories that you can think of as "classes". Each folder in src/ is associated with a different class. Each class is composed of a structure and functions associated with that structure. We declare the structure as an incomplete type in its header file so that other classes cannot access its innards directly. See http://en.wikipedia.org/wiki/Opaque_pointer. The reason we do all this is that a new user who modifies the code is much less likely to introduce subtle bugs because they can only really modify data using pre-defined functions.

* You will need to compile hdf5 with parallel support. Using macports won't give you this. I did it by simply compiling the hdf5 source using mpicc. make sure you edit the Disco Makefile according to the location of your installation.

* Viscosity is turned off when you set it to a negative number

* let N_p vary by setting NP_CONST to negative number

* turn of damping by setting DAMP_TIME to negative number

* set RLogScale and ZLogScale to large numbers to get uniform grid spacing.

* running with no inner boundary condition is currently tricky. You want to set the flag NoInnerBC=1 in your .par file, but you still need to set R_Min to be greater than 0 and big enough so that the inner face of the 1 inner ghostzone does not accidentally become negative. We will definitely figure out something better than this soon.

* Right now you need to run on a number of procs that allows you to divide the grid evenly in both r and z. This is kind of bad because when we keep the aspect ratio constant the outer procs are going to have more phi zones and our load balancing will be bad. Also, we could reduce inter-processor surface area and thus the number of ghost zones if we took NumR and NumZ into account when we figure out how to distrubute the procs. 
