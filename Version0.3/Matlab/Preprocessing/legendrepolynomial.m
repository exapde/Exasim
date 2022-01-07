function [f,fx] = legendrepolynomial(x,p)
%KOORNWINDER1D Vandermonde matrix for Legenedre polynomials in [0,1]
%   [F,FX]=KOORNWINDER(X,P)
%
%      X:         Coordinates of the points wherethe polynomials 
%                 are to be evaluated (NPOINTS)
%      P:         Maximum order of the polynomials consider. That
%                 is all polynomials of degree up to P, NPOLY=P+1
%      F:         Vandermonde matrix (NPOINTS,NPOLY)
%      FX:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. x (NPOINTS,NPOLY)
%

x = 2*x-1;

%pp=p+1;
f=zeros(numel(x),p+1);
fx=zeros(numel(x),p+1);

for ii=0:p
    pp = jacobi(ii,0,0);
% Normalization factor to ensure integration to one
    %pp = pp*sqrt(2*ii+1);  
    dpp = polyder(pp); 
    pval  = polyval(pp,x);
    dpval = polyval(dpp,x);
    f(:,ii+1) = pval;
    fx(:,ii+1) = dpval;
end

fx = 2*fx;
