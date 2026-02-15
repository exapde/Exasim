cdir = pwd();

cases{1} = cdir + "/Poisson/poisson1d/";
cases{2} = cdir + "/Poisson/poisson2d/";
cases{3} = cdir + "/Poisson/poisson3d/";
cases{4} = cdir + "/Poisson/coupledproblem/";
cases{5} = cdir + "/Poisson/periodic/";
cases{6} = cdir + "/Poisson/Axisymmetric/";
cases{7} = cdir + "/Poisson/Lshape/";
cases{8} = cdir + "/Poisson/EquationOfStates/";
cases{9} = cdir + "/Poisson/Nonlinear/";
cases{10} = cdir + "/Poisson/nonlinearheatflux/";
cases{11} = cdir + "/Poisson/Orion/";
cases{12} = cdir + "/Poisson/poisson_streamer/";
cases{13} = cdir + "/Poisson/Cone/";
cases{14} = cdir + "/Burgers/Burgers2d/";
cases{15} = cdir + "/Advection/GaussianRotating/";
cases{16} = cdir + "/ConvectionDiffusion/1D/";
cases{17} = cdir + "/ConvectionDiffusion/BoundaryLayer/";
cases{18} = cdir + "/ConvectionDiffusion/system/";
cases{19} = cdir + "/Euler/naca0012/";
cases{20} = cdir + "/Euler/EulerVortex/";
cases{21} = cdir + "/GSI/heat1d/";
cases{22} = cdir + "/GSI/concentration1d/";
cases{23} = cdir + "/GSI/concentration1d_Hf/";
cases{24} = cdir + "/GSI/heatconcentration1d/";
cases{25} = cdir + "/GSI/heatconcentration1d_HfHfO2/";
cases{26} = cdir + "/HeatEquation/warmup/";
cases{27} = cdir + "/HeatEquation/cooldown/";
cases{28} = cdir + "/LinearElasticity/beam2d/";
cases{29} = cdir + "/LinearElasticity/beam3d/";
cases{30} = cdir + "/LinearElasticity/cookmembrane/";
cases{31} = cdir + "/MHD/MagneticVortex/";
cases{32} = cdir + "/MongeAmpere/square/";
cases{33} = cdir + "/MongeAmpere/square2d/";
cases{34} = cdir + "/NavierStokes/naca2d/";
cases{35} = cdir + "/NavierStokes/flatplate2d/";
cases{36} = cdir + "/NavierStokes/hypersoniccylinder_mach8/";
cases{37} = cdir + "/NavierStokes/flaredplate2d/";
cases{38} = cdir + "/NavierStokes/nshtmach8/";
cases{39} = cdir + "/NavierStokes/naca2d_unsteady/";
cases{40} = cdir + "/NavierStokes/TaylorGreenVortex/";
cases{41} = cdir + "/NonlinearElasticity/compressionblock/";
cases{42} = cdir + "/NonlinearElasticity/cookmembrane/";
cases{43} = cdir + "/NonlinearElasticity/strechingblock/";
cases{44} = cdir + "/NonlinearElasticity/cubeblock/";
cases{45} = cdir + "/RANS/naca0012/";
cases{46} = cdir + "/RANS/flatplate/";
cases{47} = cdir + "/WaveEquation/Membrane/";
cases{48} = cdir + "/WaveEquation/Scattering/";
cases{49} = cdir + "/Stokes/Kovasznay/";
cases{50} = cdir + "/Stokes/square/";
cases{51} = cdir + "/ShallowWater/BickleyJet/";

unittests = [1:12 14:39 41:51];

for jjj = 1:length(unittests)
  iii = unittests(jjj);
  cd(cases{iii});
  clearvars -except cases unittests iii jjj
  disp("RUNNING " + cases{iii});
  run(cases{iii} + "pdeapp.m");
end


