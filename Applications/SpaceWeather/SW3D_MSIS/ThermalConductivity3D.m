function kappa = ThermalConductivity3D(x, mu, eta)
    
%Parameters that need to be defined
nspecies = 4;
nWaves = 37;

%Position
r = sqrt(x(1)^2 + x(2)^2 + x(3)^2); 

R0 = mu(16);
H0 = mu(18);

%Compute weighted density compositions (n_i/rho = Chi/mi)
Chi = sym(zeros(nspecies,1));
Chi(1) = 1.0;
for iSpecies = 2:nspecies
    coeffsDensity = eta((3+nspecies)*nWaves+4*(iSpecies-2)+1:(3+nspecies)*nWaves+4*(iSpecies-1));
    Chi(iSpecies) = coeffsDensity(1)*exp(coeffsDensity(2)*(r-R0)*H0) + coeffsDensity(3)*exp(coeffsDensity(4)*(r-R0)*H0);
    Chi(1) = Chi(1) - Chi(iSpecies);
end

%Compute thermal conductivity
kappa = sym(0.0);
ckappai = eta((3+nspecies)*nWaves+4*(nspecies-1)+nspecies+1:(3+nspecies)*nWaves+4*(nspecies-1)+2*nspecies);

for iSpecies = 1:nspecies
    kappa = kappa + ckappai(iSpecies)*Chi(iSpecies);
end

