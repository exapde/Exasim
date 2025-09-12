function [kcute,waveno,randno,lemax] = stginit(hx, hyz, dwall, lt, nu, eps, alpha)

% spectral maximum wavenumber and cut-off wavenumber
[kcute, lemax] = calculatek(hx, hyz, dwall, lt, nu, eps);

% wavenumber distributions
beta = 0.5;
kemin = min(kcute(:,1));
kcutmax = max(kcute(:,2));
waveno = wavenumbers(alpha, beta, kemin, kcutmax);

% random numbers for N modes
N = size(waveno,1);
randno = randomgen(N);





