function [stgdata,Ek,Em,dwaveno,waveno,pars] = stghomogeneousturbulence(minLength, turbLengthFactor, visc, turbIntensity, U, N)

[fluctamp,pars,Ek,Em,awaveno,dwaveno,waveno] = homogeneousturbulence(minLength, turbLengthFactor*minLength, visc, turbIntensity*U, N);
fluctamp = fluctamp/U;

randno = randomgen(length(awaveno));
randno(:,8) = 0;
for n=1:length(awaveno)
    c = turbIntensity*awaveno(n);
    %randno(n,8) = normrnd(c,c);
    randno(n,8) = c + c * randn();
end

stgdata = [awaveno fluctamp randno];
