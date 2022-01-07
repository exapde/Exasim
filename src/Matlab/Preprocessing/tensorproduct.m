function [f,fx,fy,fz] = tensorproduct(x,porder)
%TENSORPRODUCT Vandermonde matrices for tensor product polynomials in [0,1]^d
%
%   [F,FX,FY,FZ] = TENSORPRODUCT(X,P)
%
%      X:         Coordinates of the points wherethe polynomials 
%                 are to be evaluated (npoints,dim)
%      PORDER:    Maximum order of the polynomials consider. That
%                 is all polynomials of complete degree up to p,
%                 npoly = (PORDER+1)*(PORDER+1)*(PORDER+1)
%      F:         Vandermonde matrix (npoints,npoly)
%      FX:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. x (npoints,npoly)
%      FY:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. y (npoints,npoly)
%      FZ:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. z (npoints,npoly)
%

[n,dim] = size(x);

switch (dim)
    case 1 % 1D 
        [f,fx] = koornwinder(x,porder);  % Legendre basis
        fy = [];
        fz = [];
    case 2 % 2D
        if length(porder)==1
            porder = [porder porder];
        end
        [g1,gx] = koornwinder(x(:,1),porder(1)); % Legendre basis in x direction     
        [g2,gy] = koornwinder(x(:,2),porder(2)); % Legendre basis in y direction         
        f  = zeros(n,prod(porder+1));
        fx = f;
        fy = f;        
        fz = [];
        % perform tensor product to obtain the shape functions and their 
        % derivatives on the unit square
        for ii=1:n
            f(ii,:) =  kron(g2(ii,:),g1(ii,:));
            fx(ii,:) = kron(g2(ii,:),gx(ii,:));
            fy(ii,:) = kron(gy(ii,:),g1(ii,:));
        end
    case 3
        if length(porder)==1
            porder=[porder porder porder];
        end
        [g1,gx]=koornwinder(x(:,1),porder(1)); % Legendre basis in x direction         
        [g2,gy]=koornwinder(x(:,2),porder(2)); % Legendre basis in y direction             
        [g3,gz]=koornwinder(x(:,3),porder(3)); % Legendre basis in z direction             
        f  = zeros(n,prod(porder+1));
        fx = f;
        fy = f;
        fz = f;
        % perform tensor product to obtain the shape functions and their 
        % derivatives on the unit cube
        for ii=1:n
            f(ii,:) =  kron(g3(ii,:),kron(g2(ii,:),g1(ii,:)));
            fx(ii,:) = kron(g3(ii,:),kron(g2(ii,:),gx(ii,:)));
            fy(ii,:) = kron(g3(ii,:),kron(gy(ii,:),g1(ii,:)));
            fz(ii,:) = kron(gz(ii,:),kron(g2(ii,:),g1(ii,:)));
        end        
    otherwise
        error('Only can handle dim=1, dim=2 or dim=3');
end
