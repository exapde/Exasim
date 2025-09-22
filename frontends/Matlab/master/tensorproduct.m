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
%         for ii=1:n
%             f(ii,:) =  kron(g2(ii,:),g1(ii,:));
%             fx(ii,:) = kron(g2(ii,:),gx(ii,:));
%             fy(ii,:) = kron(gy(ii,:),g1(ii,:));
%         end
        m1 = size(g1,2);
        m2 = size(g2,2);        
        for i = 1:m1
          for j= 1:m2
            f(:, (i-1)*m2 + j) = g2(:,i).*g1(:,j);
            fx(:, (i-1)*m2 + j) = g2(:,i).*gx(:,j);
            fy(:, (i-1)*m2 + j) = gy(:,i).*g1(:,j);
          end
        end        
%         ft  = zeros(n,prod(porder+1));
%         ftx = f;
%         fty = f;        
%         m1 = size(g1,2);
%         m2 = size(g2,2);
%         for i = 1:m1
%           for j= 1:m2
%             ft(:, (i-1)*m2 + j) = g2(:,i).*g1(:,j);
%             ftx(:, (i-1)*m2 + j) = g2(:,i).*gx(:,j);
%             fty(:, (i-1)*m2 + j) = gy(:,i).*g1(:,j);
%           end
%         end
%         max(abs(f(:)-ft(:)))
%         max(abs(fx(:)-ftx(:)))
%         max(abs(fy(:)-fty(:)))        
%         pause
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
