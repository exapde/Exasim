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

for i = 1:length(cases)
  cd(cases{i})
  disp("RUNNING " + cases{i});
  run(cases{i} + "pdeapp.m");
end
