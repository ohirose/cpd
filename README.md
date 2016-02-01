
# Coherent Point Drift
This software is a simple implementation of the point-set registration algorithm known as Coherent   
Point Drift (CPD) invented by Andriy Myronenko and Xubo Song (2010). All of the source codes were    
written by ohirose in C language except for functions for solving linear equations, singular value   
decomposition, and eigendecomposition. Algorithm details are available in their article "Point Set  
Registration: Coherent Point Drift, IEEE TPAMI, 32(12), 2262--2275, 2010.    

## INSTALL
**Linux**: Launch terminal, move to the source file directory, and type `make`.   
**Mac**: Install Xcodes and command-line tools. After this, do the same procedure as that for Linux.   
**Windows**:  

1. Install [MinGW](http://www.mingw.org/) (+MSYS).  
2. Download [here](http://icl.cs.utk.edu/lapack-for-windows/lapack/)
   the following pre-built libraries for Windows:  
  * libgcc\_s\_dw2-1.dll 
  * libquadmath-0.dll
  * libgfortran-3.dll
  * libblas.dll
  * liblapack.lib
  * liblapack.dll
3. Move the libraries to the CPD source directory.  
4. Do the same procedure as that for Linux.  

## USAGE   
  `./cpd <mode> <X> <Y> (+ options)`

  **MODE**:   
  At least one of characters r, a, and c must be included in `mode`. Optionally, `m` and `q` which    
  specify print options can be attatched.        
  `r`: rigid, `a`: affine, `c`: cpd, `q`: quiet, `m`: memorize optimization process.   

  **INPUT**:    
  `X`: point set 1, reference points.  
  `Y`: point set 2, floating points.  

  **OPTIONs**:   
  Options must be added *after* the arguments. If the parameter file is set as the argument of `-p`,    
  other parameters specified by command-line options are ignored.   
  `-n` nloop, `-w` omega, `-l` lambda, `-b` beta, `-r` rank, `-z` zscale, `-p` parameter file,   
  `-o` output file name.  
  If the option `-v` with no argument is specified, the version information is printed.

  **EXAMPLE**:   
  `./cpd raci X.txt Y.txt -w 0.5 -l 1.5 -b 0.9 -z 3.5 -n 2000`                       


