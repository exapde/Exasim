currentdir = pwd();

familyNames = {
  "Poisson"
  "Advection"
  "Burgers"
  "ConvectionDiffusion"
  "HeatEquation"
  "Stokes"
  "LinearElasticity"
  "NonlinearElasticity"
  "WaveEquation"
  "MHD"
  "Euler"
  "GSI"
  "NavierStokes"
};

familyCases = {
  {"/Poisson/Nonlinear2/", "/Poisson/Periodic/", "/Poisson/isoq/", ...
   "/Poisson/sharpb2/", "/Poisson/multipleequations/", ...
   "/Poisson/poisson2d/", "/Poisson/poisson3d/"}
  {"/Advection/GaussianTransport/"}
  {"/Burgers/Burgers2d/"}
  {"/ConvectionDiffusion/1D/", "/ConvectionDiffusion/BoundaryLayer/", ...
   "/ConvectionDiffusion/system/"}
  {"/HeatEquation/warmup/"}
  {"/Stokes/Kovasznay/"}
  {"/LinearElasticity/beam2d/", "/LinearElasticity/beam3d/"}
  {"/NonlinearElasticity/compressionblock/", "/NonlinearElasticity/strechingblock/"}
  {"/WaveEquation/Membrane/"}
  {"/MHD/MagneticVortex/"}
  {"/Euler/EulerVortex/"}
  {"/GSI/heatconcentration1d_HfHfO2/"}
  {"/NavierStokes/naca2d/"}
};

failures = {};
familyPass = zeros(length(familyNames),1);
familyFail = zeros(length(familyNames),1);
total = 0;

for iii = 1:length(familyNames)
  family = familyNames{iii};
  disp("FAMILY " + family);

  for jjj = 1:length(familyCases{iii})
    caseRel = familyCases{iii}{jjj};
    caseDir = currentdir + caseRel;
    cd(caseDir);
    clearvars -except currentdir familyNames familyCases failures familyPass familyFail total iii jjj family caseRel caseDir
    disp("RUNNING [" + family + "] " + caseDir);

    try
      run(caseDir + "pdeapp.m");
      disp("PASS [" + family + "] " + caseDir);
      familyPass(iii) = familyPass(iii) + 1;
    catch ME
      disp("FAIL [" + family + "] " + caseDir);
      disp(ME.message);
      familyFail(iii) = familyFail(iii) + 1;
      failures{end+1,1} = family;
      failures{end,2} = caseDir;
      failures{end,3} = ME.message;
    end

    total = total + 1;
  end
end

cd(currentdir);

nfail = size(failures,1);
npass = total - nfail;
disp("Quick unit-test summary: " + string(npass) + " passed, " + string(nfail) + " failed.");

for iii = 1:length(familyNames)
  nrun = familyPass(iii) + familyFail(iii);
  if nrun > 0
    disp("  " + familyNames{iii} + ": " + string(familyPass(iii)) + " passed, " + ...
         string(familyFail(iii)) + " failed.");
  end
end

if nfail > 0
  for iii = 1:nfail
    disp("FAILED [" + failures{iii,1} + "]: " + failures{iii,2});
    disp("  " + failures{iii,3});
  end
end
