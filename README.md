<p align="center">
<img src="Documentation/exasimlogosmall.png">
</p>

# Generating Discontinuous Galerkin Codes For Extreme Scalable Simulations
Exasim is an open-source software for generating discontinuous Galerkin codes to numerically solve parametrized partial differential equations (PDEs) on different computing platforms with distributed memory.  It combines high-level languages and low-level languages to easily construct parametrized PDE models and automatically produce high-performance C++ codes. The construction of parametrized PDE models and the generation of the stand-alone C++ production code are handled by high-level languages, while the production code itself can run on various machines, from laptops to the largest supercomputers, with both CPU and Nvidia GPU processors. 

What make Exasim unique are the following distinctive features:

   - Solve a wide variety of PDEs in fluid and solid mechanics, and electromagnetism
   - Generate stand-alone C++ production code via the mathematical expressions of the PDEs
   - Implement high-order DG methods including local DG and hybridized DG methods
   - Implement diagonally implicit Runge-Kutta methods 
   - Implement parallel Newton-GMRES solvers and scalable preconditioners  
   - Employ Kokkos to provide full GPU functionality for all code components from discretization schemes to iterative solvers
   - Provide interfaces to Julia, Python, and Matlab. 
   
After downloading the Exasim source code, you can run numorous examples provided in Exasim/Applications by executing pdeapp.jl in Julia, pdeapp.py in Python, or pdeapp.m in Matlab. If Exasim needs external packages, it will install them on the fly. See the documentation https://github.com/exapde/Exasim/blob/master/doc/Exasim.pdf for more details. 

# Installation 

As Exasim generates and compiles stand-alone C++ code on the fly, Exasim does not require installation. However, Exasim needs (required) C++ compiler, (required) Blas/Lapack libaries, (optional) MPI library, (optional) Gmesh for mesh generation, (optional) METIS for mesh partitioning, (optional) Paraview for visualization, and (optional) CUDA Toolkit for Nvidia GPUs. Although these external packages can be installed by running install.jl in Julia, install.py in Python, or install.m in Matlab, it is highly recommended to use Exasim without installing external packages. Exasim will install a required package on the fly if it can not find the package on your system.

# Applications

Exasim produces C++ Code to solve a wide variety of parametrized partial differential equations from first-order, second-order elliptic, parabolic, hyperbolic PDEs, to higher-order PDEs. Many examples are provided in Exasim/Applications to illustrate how to use Exasim for solving Poisson equation, wave equation, heat equation, advection, convection-diffusion, elasticity, Euler equations, Navier-Stokes equations, and MHD equations. See https://github.com/exapde/Exasim/blob/master/Applications/ShallowWater/BickleyJet/BickleyJet.pdf for simulation results of the Bickley Jet problem.


# Publications
[1] Under-Resolved Direct Numerical Simulation of Transonic Buffet Using an Implicit Discontinuous Galerkin Method. https://github.com/exapde/Exasim/blob/master/Documentation/AIAA2021paper.pdf. The numerical results and the source codes are located at https://github.com/exapde/Exasim/blob/master/Applications/NavierStokes/OAT15A.
