function p = jacobi(n,a,b)
%JACOBI returns the Jacobi polynomial of order n, 
%            P_n^{a,b}(x), for x in [-1,1]
%
%   [P]=JACOBI(N,A,B)
%
%      P:         vector of lenght n+1 whose elements are the 
%                 coefficients of the polynomial in descending powers.
%
%             y = p(1)*x^n + p(2)*x^(n-1) + ... + p(n)*x + p(n+1)
%             (y can be evaluated using the matlab function polyval)
%
%      A:         Jacobi polynomial parameter   
%      B:         Jacobi polynomial parameter
%
%    Recall that P_n^{0,0}(x) = L_n(x) are the Legenedre polynomials
%                P_n^{-1/2,-1/2}(x) = T_n(X) are the Chevychev polynomials
%
%
p0 = [1];
p1 = [(a+b+2)/2, (a-b)/2];

if n==0
    p = p0;
elseif n==1
    p = p1;
elseif n>1        
    for i=1:n-1
        a1 = 2*(i+1)*(i+a+b+1)*(2*i+a+b);
        a2 = (2*i+a+b+1)*(a*a-b*b);
        a3 = (2*i+a+b)*(2*i+a+b+1)*(2*i+a+b+2);
        a4 = 2*(i+a)*(i+b)*(2*i+a+b+2);
        p2 = conv([a3,a2],p1);
        q = zeros(1,i+2);
        q(3:i+2) = a4*p0;
        p = (p2-q)/a1;
        p0 = p1;
        p1 = p;
    end
end