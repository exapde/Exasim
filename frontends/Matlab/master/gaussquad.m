function [x,w] = gaussquad(pgauss,dim,elemtype)
%GAUSSQUAD Calculates the Gauss integration points for simplex and tensor
%          elements
%
%   [X,W]=GAUSSQUAD(PGAUSS, DIM, ELEMTYPE)
%
%      PGAUSS:    Order of the polynomial integrated exactly
%      DIM:       1D, 2D or 3D
%      ELEMTYPE:  Flag determining element type
%                 Flag = 0 simplex elements (default)
%                 Flag = 1 tensor elements
%
%      X:         Coordinates of the integration points 
%      W:         Weights  


if nargin<2, dim = 1; end

switch dim
    case 0
        x = zeros(1,0); w = 1;
    case 1 % 1D
        [x,w] = gaussquad1d(pgauss);
    case 2 % 2D
        if elemtype == 0
            [x,w] = gaussquad2d(pgauss(1));      % tri
        elseif elemtype == 1
            [x,w] = gaussquadts(pgauss,dim);  % quad
        else
            error('Element type is not implemented');
        end
    case 3 % 3D
        if elemtype == 0
            [x,w] = gaussquad3d(pgauss(1));      % tet
        elseif elemtype == 1
            [x,w] = gaussquadts(pgauss,dim);  % hex
        elseif elemtype == 2
            [x,w] = gaussquadprism(pgauss);  % prism
        elseif elemtype == 3
            [x,w] = gaussquadpyramid(pgauss);  % pyramid    
        else
            error('Element type is not implemented');    
        end
    otherwise
        error('Dimension not implemented.');
end


function [x,w] = gaussquadts(pgauss,dim)
%GAUSSQUADTS Calculates the Gauss integration points for simplex elements
%
%   [X,W]=GAUSSQUADTS(PGAUSS, DIM)
%
%      PGAUSS:    Order of the polynomial integrated exactly
%      DIM:       1D, 2D or 3D
%
%      X:         Coordinates of the integration points 
%      W:         Weights  

if length(pgauss)==1
    pgauss = pgauss*ones(dim,1);
end

% quadratures in 1D
for i = 1:dim
    [x1d{i},w1d{i}] = gaussquad1d(pgauss(i));
end

% perform tensor products to obtain quadrature rules in nD
switch dim
    case 1
        x = x1d{1};
        w = w1d{1};
    case 2
        [x2,y2] = ndgrid(x1d{1},x1d{2});
        x = [x2(:),y2(:)];        
        w = kron(w1d{2},w1d{1});
    case 3        
        [x2,y2,z2] = ndgrid(x1d{1},x1d{2},x1d{3});
        x = [x2(:),y2(:),z2(:)];        
        w = kron(kron(w1d{3},w1d{2}),w1d{1});
    otherwise
    error('Dimension not implemented.');
end

function [x,w] = gaussquadprism(pgauss)

if length(pgauss)==1
    pgauss = pgauss*ones(2,1);
end

[xtri,wtri] = gaussquad2d(pgauss(1));      % tri
[x1d,w1d] = gaussquad1d(pgauss(2));
ntri = size(xtri,1);
n1d = length(x1d);
x = zeros(ntri*n1d,3);
for i=1:n1d
    ind = (1:1:ntri) + (i-1)*ntri;
    x(ind,1:2) = xtri;
    x(ind,3) = x1d(i);
end
w = kron(w1d,wtri);

