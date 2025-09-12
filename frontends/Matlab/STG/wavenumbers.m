function [waveno, kmin, kmax] = wavenumbers(alpha, beta, kemin, kcutmax)

% minimum wavenumber
kmin = beta*kemin; 

% maximum wavenumber
kmax = 1.5*kcutmax;

% number of modes
N = 1 + ceil( log(kmax / kmin) / log(1.0 + alpha));

% calculate wavenumber distribution
waveno = zeros(N,1);
waveno(1) = kmin;
for n = 2:N
    waveno(n) = waveno(n-1)*(1+alpha);
end

% average wavenumber distribution
awaveno = 0.5*(waveno(2:end) + waveno(1:end-1));

% difference wavenumber distribution
dwaveno = waveno(2:end) - waveno(1:end-1);

% average and difference
waveno = [awaveno dwaveno];







