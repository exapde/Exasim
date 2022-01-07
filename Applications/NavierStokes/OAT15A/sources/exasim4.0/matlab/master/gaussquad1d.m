function [x,w]=gaussquad1d(pgauss)
%GAUSSQUAD1D Calculates the Gauss integration points in 1D for [0,1]
%   [X,W]=GAUSSQUAD1D(PGAUSS)
%
%      X:         Coordinates of the integration points 
%      W:         Weights  
%      PGAUSS:    Order of the polynomial integrated exactly
%

% % Our original implementation
% n=ceil((pgauss+1)/2);
% P=jacobi(n,0,0);
% x=sort(roots(P));
% 
% A=zeros(n,n);
% for i=1:n
%   P = jacobi(i-1,0,0);
%   A(i,:)=polyval(P,x)';
% end
% w=A\[2;zeros(n-1,1)];
% x=(x+1)/2;
% w=w/2;

% Trefethen implementation:
n=ceil((pgauss-1)/2);
beta = .5./sqrt(1-(2*(1:n)).^(-2)); % 3-term recurrence coeffs
T = diag(beta,1) + diag(beta,-1); % Jacobi matrix
[V,D] = eig(T); % eigenvalue decomposition
x = diag(D); [x,i] = sort(x); % nodes (= Legendre points)
w = 2*V(1,i).^2; % weights
% normalize to [0 1]
w = w'/2;
x = .5*(x + 1);
