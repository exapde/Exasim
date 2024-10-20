function XT=block2tets(X,ind)
%BLOCK2TETS  Subdivides a block into tets - allows for a parity index.
%
%   XT=block2tets(X,ind)

if nargin<2, ind=1; end

dim = size(X,1);

switch dim
    case 1
        XT = [X(1:end-1);X(2:end)];
        XT = reshape(XT,1,2,[]);
    case 2
        nx = size(X,2); ny = size(X,3);
        XT = zeros(2,3,2*(nx-1)*(ny-1));
        in = ind;
        ik = 1;
        for j=1:ny-1
            iy = in;
            for i=1:nx-1
                XL = Gtris( X(:,i,j), X(:,i+1,j), X(:,i,j+1), X(:,i+1,j+1), in);
                XT(1:2,1:3,ik:ik+1) = XL;
                ik = ik+2;
                in = -in;
            end
            in = -iy;
        end
    case 3
        nx = size(X,2); ny = size(X,3); nz = size(X,4);
        XT = zeros(3,4,5*(nx-1)*(ny-1)*(nz-1));
        in = ind;
        ik = 1;
        for k=1:nz-1
            iz = in;
            for j=1:ny-1
                iy = in;
                for i=1:nx-1
                    XL = Gtets( X(:,i,j,k  ), X(:,i+1,j,k  ), X(:,i,j+1,k  ), X(:,i+1,j+1,k  ), ...
                                X(:,i,j,k+1), X(:,i+1,j,k+1), X(:,i,j+1,k+1), X(:,i+1,j+1,k+1), in);
                    XT(1:3,1:4,ik:ik+4) = XL;
                    ik = ik+5;
                    in = -in;
                end
                in = -iy;
            end
            in = -iz;
        end
end

function XL = Gtris( P0, P1, P2, P3, ind)

XL = zeros(2,3,2);

if ind == 1
    XL(:,1,1)= P0; XL(:,2,1)= P1; XL(:,3,1)= P3;
    XL(:,1,2)= P0; XL(:,2,2)= P3; XL(:,3,2)= P2;
elseif ind == -1
    XL(:,1,1)= P0; XL(:,2,1)= P1; XL(:,3,1)= P2;
    XL(:,1,2)= P1; XL(:,2,2)= P3; XL(:,3,2)= P2;
end


function XL = Gtets( P0, P1, P2, P3, P4, P5, P6, P7, ind)

XL = zeros(3,4,5);

if ind == 1
    XL(:,1,1)= P0; XL(:,2,1)= P1; XL(:,3,1)= P3; XL(:,4,1)= P5;
    XL(:,1,2)= P0; XL(:,2,2)= P3; XL(:,3,2)= P2; XL(:,4,2)= P6;
    XL(:,1,3)= P4; XL(:,2,3)= P6; XL(:,3,3)= P5; XL(:,4,3)= P0;
    XL(:,1,4)= P7; XL(:,2,4)= P5; XL(:,3,4)= P6; XL(:,4,4)= P3;
    XL(:,1,5)= P0; XL(:,2,5)= P5; XL(:,3,5)= P3; XL(:,4,5)= P6;    
elseif ind == -1
    XL(:,1,1)= P0; XL(:,2,1)= P1; XL(:,3,1)= P2; XL(:,4,1)= P4;
    XL(:,1,2)= P1; XL(:,2,2)= P3; XL(:,3,2)= P2; XL(:,4,2)= P7;
    XL(:,1,3)= P4; XL(:,2,3)= P7; XL(:,3,3)= P5; XL(:,4,3)= P1;
    XL(:,1,4)= P4; XL(:,2,4)= P6; XL(:,3,4)= P7; XL(:,4,4)= P2;
    XL(:,1,5)= P1; XL(:,2,5)= P4; XL(:,3,5)= P7; XL(:,4,5)= P2;  
end
