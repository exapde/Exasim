function [waveno, awaveno, dwaveno] = wavenumberdistr(N, kmin, kmax)

% wave number ratio
rat = exp(log((kmax/kmin))/(N-1));

% calculate wavenumber distribution
waveno = zeros(N,1);
waveno(1) = kmin;
for n = 1:N-1
    waveno(n+1) = waveno(n)*rat;
end
% deltak = (log(kmax) - log(kmin))/(N-1);
% waveno = exp(log(kmin) + (0:(N-1))'*(deltak));
%waveno = linspace(kmin,kmax,N);

% average wavenumber distribution
awaveno = 0.5*(waveno(2:end) + waveno(1:end-1));

% difference wavenumber distribution
dwaveno = waveno(2:end) - waveno(1:end-1);








