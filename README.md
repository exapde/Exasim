<p align="center">
<img src="doc/exasimlogosmall.png">
</p>

# Generating Discontinuous Galerkin Codes For Extreme Scalable Simulations
Exasim is an open-source software for generating high-order discontinuous Galerkin (DG) codes to numerically solve parametrized partial differential equations (PDEs) on different computing platforms with distributed memory.  It combines high-level languages and low-level languages to easily construct parametrized PDE models and automatically produce high-performance C++ codes. The construction of parametrized PDE models and the generation of the stand-alone C++ production code are handled by high-level languages, while the production code itself can run on various machines, from laptops to the largest supercomputers, with  AMD and Nvidia GPU processors. Exasim has the following capabilities:

   - Solve a wide variety of PDEs in fluid mechanics, solid mechanics, electromagnetism, and multi-physics models, in 1D, 2D, and 3D
   - Generate stand-alone C++ production code via the mathematical expressions of the PDEs
   - Implement local DG and hybridized DG methods for spatial discretization
   - Implement diagonally implicit Runge-Kutta methods for temporal discretization
   - Implement parallel Newton-GMRES solvers and scalable preconditioners using reduced basis method and polynomial preconditioners
   - Implement monolithic multi-physics solvers for the HDG discretization    
   - Employ Kokkos to provide full GPU functionality for all code components from discretization schemes to iterative solvers
   - Provide auto-gen tools to calculate thermodynamic, transport, chemistry, and energy transfer properties for chemically-reacting flows 
   - Provide application interfaces to Julia, Python, and Matlab. 
   
After downloading the source code, please make sure that the name of the folder is `Exasim`. If it has a different name, please rename it to `Exasim`. Please make sure that the directory containing the folder Exasim does not have any white space, because Kokkos libraries can not be compiled properly in such case. See [the documentation](https://github.com/exapde/Exasim/blob/master/doc/Exasim.pdf) for more details. 

To deploy, compile, and run Exasim on **HPC systems**, please follow the intructions in [the hpc manual](https://github.com/exapde/Exasim/blob/master/install/hpc.txt).

## Installation 

Exasim needs Kokkos (required), Blas/Lapack libaries (required), MPI library (required), Gmesh for mesh generation (optional), METIS for mesh partitioning (optional), Paraview for visualization (optional), and CUDA Toolkit (optional) to run on Nvidia GPUs. 

### Build Kokkos Libaries
Since Exasim uses Kokkos to target various computing platforms, you must build Kokkos libraries before using Exasim. To build Kokkos serial library for CPU platform:
```
cd Exasim/kokkos   
make -f Makefile.builds serial 
```

To build Kokkos CUDA library for Nvidia GPU platform:
```
cd Exasim/kokkos
make -f Makefile.builds cuda   
```
To build Kokkos HIP library for AMD GPU platform:
```
cd Exasim/kokkos
make -f Makefile.builds hip   
```

### Build Metis and ParMetis Libraries

To build Metis and ParMetis libaries:
```
cd Exasim/metis
make metis
```

### Build Text2Code Executable

To use Text2Code as a code generator and preprocessor in Exasim:
```
cd Exasim/text2code
make text2code
``` 

Text2Code produces faster, cleaner code and removes the need for MATLAB/Julia/Python at runtime. 

### Build Exasim Executables

After installing Text2Code successfully, please procceed installing Exasim as follows
```
cd Exasim/build
cmake -D EXASIM_NOMPI=ON -D EXASIM_MPI=ON -D EXASIM_CUDA=OFF -D EXASIM_HIP=OFF -D WITH_TEXT2CODE=ON -D WITH_PARMETIS=ON ../install 
cmake --build .
``` 

EXASIM_CUDA=ON switches to the CUDA backend (ensure kokkos/buildcuda exists). Similarly, EXASIM_HIP=ON switches to the HIP backend. It will produce Exasim's executable programs in Exasim/build, which are **cput2cEXASIM** for CPU platform on one core, **cpumpit2cEXASIM** for CPU platform on many cores, **gput2cEXASIM** for CUDA/HIP platform on one GPU, and **gpumpit2cEXASIM** for CUDA/HIP platform on many GPUs. 

## Examples

Exasim produces C++ Code to solve a wide variety of parametrized partial differential equations from first-order, second-order elliptic, parabolic, hyperbolic PDEs, to higher-order PDEs. Many examples are provided in `Exasim/examples` to illustrate how to use Exasim for solving Poisson equation, wave equation, heat equation, advection, convection-diffusion, linear elasticity, nonlinear elasticity, Euler equations, Navier-Stokes equations, and MHD equations.   

To run any example with Julia, type the following line and hit return

```
   julia> include("pdeapp.jl")
```

To run any example with Python,  type the following line and hit return

```
   > > > exec(open("pdeapp.py").read())
```

To run any example with Matlab, type the following line and hit return

```
   > >  pdeapp
```

Exasim produces two folders in the `Exasim/build` directory. The `datain` folder contains input files which store the master, mesh and initial solution. The `dataout` folder contains the output files which store the numerical solution of the PDE model defined in the `pdeapp` script.

To run an example using Text2Code, open a terminal and perform the following commands

```
cd /path/to/Exasim/examples/<example>
/path/to/Exasim/build/text2code pdeapp.txt
/path/to/Exasim/build/cput2cEXASIM pdeapp.txt                     (if you run on one CPU core)
/path/to/Exasim/build/gput2cEXASIM pdeapp.txt                     (if you run on one GPU)
mpirun -np $N /path/to/Exasim/build/cpumpit2cEXASIM pdeapp.txt    (if you run on many CPU cores)
mpirun -np $N /path/to/Exasim/build/gpumpit2cEXASIM pdeapp.txt    (if you run on many GPUs) 
```

where N is the number of processors you specify in pdeapp.txt. Make sure to set MPI and GPU environment variables appropriately on your system. If there are examples that do not have pdeapp.txt and pdemodel.txt, they can be made by making use of pdeapp.m and pdemodel.m.

## Publications
[1] Vila-Pérez, J., Van Heyningen, R. L., Nguyen, N.-C., & Peraire, J. (2022). Exasim: Generating discontinuous Galerkin codes for numerical solutions of partial differential equations on graphics processors. SoftwareX, 20, 101212. https://doi.org/10.1016/j.softx.2022.101212

[2] Hoskin, D. S., Van Heyningen, R. L., Nguyen, N. C., Vila-Pérez, J., Harris, W. L., & Peraire, J. (2024). Discontinuous Galerkin methods for hypersonic flows. Progress in Aerospace Sciences, 146, 100999. https://doi.org/10.1016/j.paerosci.2024.100999

[3] Nguyen, N. C., Terrana, S., & Peraire, J. (2022). Large-Eddy Simulation of Transonic Buffet Using Matrix-Free Discontinuous Galerkin Method. AIAA Journal, 60(5), 3060–3077. https://doi.org/10.2514/1.j060459

[4] Nguyen, N. C., & Peraire, J. (2012). Hybridizable discontinuous Galerkin methods for partial differential equations in continuum mechanics. Journal of Computational Physics, 231(18), 5955–5988. https://doi.org/10.1016/j.jcp.2012.02.033
