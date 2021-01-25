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
   - can be called from Julia, Python, and Matlab. 

# Installation 

As Exasim generates and compiles stand-alone C++ code on the fly, Exasim does not require installation. However, Exasim needs (required) C++ compiler, (required) Blas/Lapack libaries, (optional) MPI library, (optional) Gmesh for mesh generation, (optional) METIS for mesh partitioning, (optional) Paraview for visualization, and (optional) CUDA Toolkit for Nvidia GPUs. To install these packages, go to the directory Exasim/Installation and run install.jl in Julia, install.py in Python, or install.m in Matlab. See the documentation https://github.com/exapde/Exasim/tree/master/Documentation/Exasim0.3.pdf for more details. 

# Applications

Exasim produces C++ Code to solve a wide variety of parametrized partial differential equations from first-order, second-order elliptic, parabolic, hyperbolic PDEs, to higher-order PDEs. Many examples are provided to illustrate how to use Exasim for solving Poisson equation, wave equation, heat equation, advection, convection-diffusion, elasticity, Euler equations, Navier-Stokes equations, and MHD equations. 

To try out any of the provided examples, please download the Exasim source code and go to any folder under Exasim/Applications and run pdeapp.jl in Julia, pdeapp.py in Python, or pdeapp.m in Matlab. See https://github.com/exapde/Exasim/blob/master/Applications/ShallowWater/BickleyJet/BickleyJet.pdf for simulation results of the Bickley Jet problem.

# Publications
[1] A Matrix-Free Implicit Discontinuous Galerkin Method for Large Eddy Simulation of Transonic Buffet at High Reynolds Number on Graphics Processors. https://github.com/exapde/Exasim/tree/master/Documentation/AIAA_Journal_2020_OAT15A.pdf 
