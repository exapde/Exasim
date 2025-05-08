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
   - Provide application programming interfaces to Julia, Python, and Matlab. 
   
After downloading the source code, please make sure that the name of the folder is `Exasim`. If it has a different name, please rename it to `Exasim`. Please make sure that the directory containing the folder Exasim does not have any white space, because Kokkos libraries can not be compiled properly in such case. See [the documentation](https://github.com/exapde/Exasim/blob/master/doc/Exasim.pdf) for more details. 

To deploy, compile, and run Exasim on **HPC systems**, please follow the intructions in [the hpc manual](https://github.com/exapde/Exasim/blob/master/install/hpc.txt).

# Installation 

Exasim needs Kokkos (required), Blas/Lapack libaries (required), MPI library (required), Gmesh for mesh generation (optional), METIS for mesh partitioning (optional), Paraview for visualization (optional), and CUDA Toolkit (optional) to run on Nvidia GPUs. These external packages can be installed by running install.jl in Julia, install.py in Python, or install.m in Matlab.

As Exasim generates and compiles stand-alone C++ code on the fly, Exasim does not require installation. However, since Exasim uses Kokkos to target various computing platforms, you must build Kokkos libraries before using Exasim. To build Kokkos serial library for CPU platform, please follow the below steps

```
  $ cd Exasim/kokkos   
  $ mkdir buildserial
  $ cd buildserial
  $ cmake .. -DCMAKE_INSTALL_PREFIX=../buildserial
  $ make install   
```

To build Kokkos CUDA library for Nvidia GPU platform, please follow the below steps
```
  $ cd Exasim/kokkos
  $ mkdir buildcuda
  $ cd buildcuda
  $ cmake .. -DCMAKE_CXX_COMPILER=clang++ -DKokkos_ENABLE_CUDA=ON -DCMAKE_INSTALL_PREFIX=../buildcuda
  $ make install   
```
To build Kokkos HIP library for AMD GPU platform, please follow the below steps
```
  $ cd Exasim/kokkos
  $ mkdir buildhip
  $ cd buildhip
  $ cmake .. -DCMAKE_CXX_COMPILER=hipcc -DKokkos_ENABLE_HIP=ON -DKokkos_ENABLE_ROCM=ON -DCMAKE_INSTALL_PREFIX=../buildhip
  $ make install   
```

Once Kokkos libraries are successfully built, you can start using Exasim. To try out any of the provided examples, please go to any folder in the directory  Exasim/examples and run pdeapp.jl in Julia, pdeapp.py in Python, or pdeapp.m in Matlab. 

# Examples

Exasim produces C++ Code to solve a wide variety of parametrized partial differential equations from first-order, second-order elliptic, parabolic, hyperbolic PDEs, to higher-order PDEs. Many examples are provided in `Exasim/examples` to illustrate how to use Exasim for solving Poisson equation, wave equation, heat equation, advection, convection-diffusion, Euler equations, Navier-Stokes equations, and MHD equations. See [the Bickley Jet example](https://github.com/exapde/Exasim/blob/master/examples/ShallowWater/BickleyJet/BickleyJet.pdf) for simulation results.

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

If successful, Exasim produces an executable application and three new folders in the `build` folder. The `backend/Model` folder contains the  Kokkos source code generated by Exasim, the `build/datain` folder contains input files for the executable application, and the `build/dataout` folder contains the output files produced by running the executable application,  which stores the numerical solution of the PDE model defined in the `pdeapp` script. The name of the executable application is **cpuEXASIM** for CPU platform on one processor, **cpumpiEXASIM** for CPU platform on many processors, **gpuEXASIM** for CUDA platform on one GPU, and **gpumpiEXASIM** for CUDA platform on many GPUs. 

# Publications
[1] Vila-Pérez, J., Van Heyningen, R. L., Nguyen, N.-C., & Peraire, J. (2022). Exasim: Generating discontinuous Galerkin codes for numerical solutions of partial differential equations on graphics processors. SoftwareX, 20, 101212. https://doi.org/10.1016/j.softx.2022.101212

[2] Hoskin, D. S., Van Heyningen, R. L., Nguyen, N. C., Vila-Pérez, J., Harris, W. L., & Peraire, J. (2024). Discontinuous Galerkin methods for hypersonic flows. Progress in Aerospace Sciences, 146, 100999. https://doi.org/10.1016/j.paerosci.2024.100999

[3] Nguyen, N. C., Terrana, S., & Peraire, J. (2022). Large-Eddy Simulation of Transonic Buffet Using Matrix-Free Discontinuous Galerkin Method. AIAA Journal, 60(5), 3060–3077. https://doi.org/10.2514/1.j060459

[4] Nguyen, N. C., & Peraire, J. (2012). Hybridizable discontinuous Galerkin methods for partial differential equations in continuum mechanics. Journal of Computational Physics, 231(18), 5955–5988. https://doi.org/10.1016/j.jcp.2012.02.033
