<p align="center">
<img src="Documentation/exasimlogosmall.png">
</p>

# Building Partial Differential Equation Applications For Extreme Scalable Simulations
Exasim is an open-source software for generating discontinuous Galerkin codes and building executable applications to numerically solve  parametrized partial differential equations (PDEs) on different computing platforms with distributed memory.  It combines high-level languages  and low-level languages to easily create parametrized PDE models and automatically produce high-performance applications. The construction of parametrized PDE models and the generation of the stand-alone C++ production code are handled by high-level languages, while the production code itself can run on various machines, from laptops to the largest supercomputers, with both CPU and Nvidia GPU processors. 

# Dependency  

Exasim automatically generates and compiles stand-alone C++ code on the fly. To do that, Exasim requires a C++ compiler and Blas/Lapack libraries for building serial applications. An MPI library is required to build parallel applications. CUDA Tookit is required to build applications on Nvidia GPUs. Gmesh is used for mesh generation. METIS is needed for mesh partitioning. And Paraview is used for visualization. To install these packages, please go to the directory Exasim/Installation and run install.jl in Julia, install.py in Python, or install.m in Matlab. See the documentation https://github.com/exapde/Exasim/tree/master/Documentation/Exasim0.1.pdf for more details. 

# Language Support

Exasim is available in Julia, Python, and Matlab. 

# Partial Differential Equations

Exasim produces C++ Code to solve a wide variety of parametrized partial differential equations from first-order PDEs, second-order elliptic, parabolic, hyperbolic PDEs, to higher-order PDEs.

# Applications

Many examples are provided to illustrate how to build executable applications for solving a wide variety of PDEs including Poisson equation, wave equation, heat equation, advection, convection-diffusion, elasticity, Euler equations, Navier-Stokes equations, and MHD equations. To try out any of the provided examples, please go to any folder under Exasim/Applications and run pdeapp.jl in Julia, pdeapp.py in Python, or pdeapp.m in Matlab. 
