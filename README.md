<p align="center">
<img src="Documentation/exasimlogosmall.png">
</p>

# Generating Discontinuous Galerkin Codes For Extreme Scalable Simulations
Exasim is an open-source software for generating discontinuous Galerkin codes to numerically solve parametrized partial differential equations (PDEs) on different computing platforms with distributed memory.  It combines high-level languages and low-level languages to easily construct parametrized PDE models and automatically produce high-performance C++ codes. The construction of parametrized PDE models and the generation of the stand-alone C++ production code are handled by high-level languages, while the production code itself can run on various machines, from laptops to the largest supercomputers, with both CPU and Nvidia GPU processors. 

What make Exasim unique are the following distinctive features:

   - simplifies numerical modeling and simulation with simple scripts
   - produces implicit high-order discontinuous Galerkin solutions of a wide variety of PDEs   
   - generates stand-alone C++ production code tailored to specific applications on different platforms 
   - provides full GPU functionality, meaning that all code components from discretization schemes to iterative solvers run on GPUs.   

# Dependency  

Exasim automatically generates and compiles stand-alone C++ code on the fly. To do that, Exasim requires a C++ compiler and Blas/Lapack libraries for generating serial code. An MPI library is required to compile and run parallel code. CUDA Tookit is required to run CUDA code on Nvidia GPUs. Gmesh is used for mesh generation. METIS is needed for mesh partitioning. And Paraview is used for visualization. 

# Installation  

Exasim itself does not require installation. However, Exasim needs (required) C++ compiler, (required) Blas/Lapack libaries, (optional) MPI library, (optional) Gmesh, (optional) Paraview, and (optional) CUDA Toolkit. To install any of these packages, please go to the directory Exasim/Installation and run install.jl in Julia, install.py in Python, or install.m in Matlab. See the documentation https://github.com/exapde/Exasim/tree/master/Documentation/Exasim0.3.pdf for more details. 

# Language Support

Exasim is available in Julia, Python, and Matlab. 

# Partial Differential Equations

Exasim produces C++ Code to solve a wide variety of parametrized partial differential equations from first-order, second-order elliptic, parabolic, hyperbolic PDEs, to higher-order PDEs.

# Applications

Many examples are provided to illustrate how to build stand-alone C++ codes for solving a wide variety of PDEs including Poisson equation, wave equation, heat equation, advection, convection-diffusion, elasticity, Euler equations, Navier-Stokes equations, and MHD equations. To try out any of the provided examples, please go to any folder under Exasim/Applications and run pdeapp.jl in Julia, pdeapp.py in Python, or pdeapp.m in Matlab. See https://github.com/exapde/Exasim/blob/master/Applications/ShallowWater/BickleyJet/BickleyJet.pdf for simulation results of the Bickley Jet problem.

# Publications
[1] A Matrix-Free Implicit Discontinuous Galerkin Method for Large Eddy Simulation of Transonic Buffet at High Reynolds Number on Graphics Processors. https://github.com/exapde/Exasim/tree/master/Documentation/AIAA_Journal_2020_OAT15A.pdf 
