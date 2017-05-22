# README 

`Bathe` is a code made for studying the physical behavior of solids. It uses the finite element method to solve the
non-linear system of equations. Several schemes like Newton-Raphson or Arc-Length of load control are implemented and
they should be used depending of the problem.  

## Instalation

The main libraries used by this code are PETSc and SLEPc. Up to now only PETSc is used to solve the linear system of
equations. SLEPc will be used in the future in case of studing eigensystems.

###`BLAS & LAPACK` libraries

```bash
apt-get install libblas-dev liblapack-dev
```  

###`PETSC` library

Download it from [www.mcs.anl.gov/petsc](www.mcs.anl.gov/petsc) and do:

```bash    
   export PETSC_DIR=/path-to-petsc-installation-directory
   export PETSC_ARCH=arch-linux2-c-opt    
   LIBRARY_PATH=$LD_LIBRARY_PATH:/path-to-petsc-directory
```  
we recommend to put them in *.bashrc* or in some start up file

```bash
   tar -xvzf petsc-{version}
   ./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-fblaslapack --download-mpich
   make all 
   make test
```

###`SLEPc` library:

Download it from [http://slepc.upv.es](http://slepc.upv.es) and do:

```bash   
   export SLEPC_DIR=/path-to-slepc-installation-directory 
```

we recommend to put them in `.bashrc` or in some start up file

```bash
   tar -xvzf slepc-{version} 
   ./configure
   make all
   make test
```
###`Bathe`:

```bash
    make
```

##Running `Bathe` on multiple processors

```bash
   mpirun -np 2  bathe <inputfile.bathe>
```

##Debbuging `Bathe` on multiple processors

```bash
   mpirun -np 2  xterm -e gdb --args bathe <inputfile.bathe>
```   
-ksp_atol <abstol>  - Sets abstol

-ksp_rtol <rtol>    - Sets rtol

-ksp_divtol <dtol>  - Sets dtol

-ksp_max_it <maxits>    - Sets maxits 

-ksp_monitor 200


## Do a PDF of this  

For reading this text in a `pdf` format do:

```bash
pandoc README.md -V geometry:margin=.5in --latex-engine=xelatex -o README.pdf
```

## The future  

* Paralelization and performance evaluation 

* Scripts to test the code

* Benchmarking

* Documentation

Guido Giuntoli - [giuntoli1991@gmail.com]
