function [mw,dmdr] = weightedMass3D(x,mu,eta)

nspecies = 4;
nWaves = 37;

%Position
r = sqrt(x(1)^2 + x(2)^2 + x(3)^2);

R0 = mu(16);
H0 = mu(18);

%Compute weighted density compositions (n_i/rho = Chi/mi)
Chi = sym(zeros(nspecies,1));
dChidr = sym(zeros(nspecies,1));
Chi(1) = 1.0;
for iSpecies = 2:nspecies
    coeffsDensity = eta((3+nspecies)*nWaves+4*(iSpecies-2)+1:(3+nspecies)*nWaves+4*(iSpecies-1));
    Chi(iSpecies) = coeffsDensity(1)*exp(coeffsDensity(2)*(r-R0)*H0) + coeffsDensity(3)*exp(coeffsDensity(4)*(r-R0)*H0);
    Chi(1) = Chi(1) - Chi(iSpecies);
    
    dChidr(iSpecies) = (coeffsDensity(1)*coeffsDensity(2)*exp(coeffsDensity(2)*(r-R0)*H0) + coeffsDensity(3)*coeffsDensity(4)*exp(coeffsDensity(4)*(r-R0)*H0))*H0;
    dChidr(1) = dChidr(1) - dChidr(iSpecies);
end

mass = eta((3+nspecies)*nWaves+4*(nspecies-1)+1:(3+nspecies)*nWaves+4*(nspecies-1)+nspecies);
mw = sym(0.0);
dmdr = sym(0.0);
for iSpecies = 1:nspecies
    mw = mw + mass(iSpecies)*Chi(iSpecies);
    dmdr = dmdr + mass(iSpecies)*dChidr(iSpecies);
end

    