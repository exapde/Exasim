function [M,Minv,Smin,Smax] = meshmetric2(dgnodes,porder,elemtype,nodetype,smoothing)
%MESHMETRIC2 Get some mesh metric terms, and gradient of elements mapping
%Used for the Metric-based sensors for artificial viscosities
%
% SYNTAX: [Mi,Minv,Smin] = meshmetric2(dgnodes,porder,elemtype,nodetype)
%
% INPUTS:
%    DGNODES  - Nodes positions of the mesh
%    PORDER   - Polynomial order of the mesh
%    ELEMTYPE - Type of element (0=simplex, 1=tensor)
%    NODETYPE - Nodal arrangement type (0=regular, 1=nonregular)
%
% OUTPUTS:
%    M       - Metric Tensor of Element Jacobian (J) Mapping : M = J*J'
%    MINV    - Inverse of the Metric Tensor
%    SMIN    - Minimum Eigenvalue of the Metric Tensor
%
% SUBFUNCTIONS: voljac
%
% SEE ALSO: MESHMETRIC
%
% Author(s): NGUYEN, TERRANA
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email address: XXXXXXX@mit.edu 
% Website: http://aeroastro.mit.edu/
% November 2018

if nargin<5
    smoothing = 1;
end

% Get dimensions
[np,nd,ne] = size(dgnodes);
nm = np*ne;

% Initialization
M    = zeros(nd,nd,nm);
Minv = zeros(nd,nd,nm);
Smin = zeros(nm,1);
Smax = zeros(nm,1);

% Get derivatives at dgnodes
plocvl = mkmasternodes(porder,nd,elemtype,nodetype);
shapnl = mkshape(porder,plocvl,plocvl,elemtype);
shapnt = permute(shapnl,[2 1 3]);
[dXdXi,tm] = voljac(shapnt,dgnodes,smoothing);
if max(abs(tm(:)-dgnodes(:)))>1e-12
    error('something wrong');
end

% Create coefficients for Simplexes
A = eye(nd,nd);
if (nd==2) && (elemtype==0)
    A(1,2) = -1.0/sqrt(3.0);
    A(2,2) = 2.0/sqrt(3.0);
elseif (nd==3) && (elemtype==0)
    A(1,2) = -1.0/sqrt(3.0);
    A(2,2) = sqrt(4.0/3.0);
    A(1,3) = -sqrt(1.0/6.0);
    A(2,3) = -sqrt(1.0/6.0);
    A(3,3) = sqrt(3.0/2.0);
end

for i = 1:nm
    % Note that dXdXi is the transpose of matrix dX_i/dxi_j
    B = dXdXi(:,:,i)*A;
    M(:,:,i) = B'*B;
    Minv(:,:,i) = inv(M(:,:,i));
    %Smin(i) = eigs(squeeze(M(:,:,i)),1,'SA');
    Smin(i) = (det(dXdXi(:,:,i)))^(1/nd);
    Smax(i) = sqrt(eigs(squeeze(M(:,:,i)),1,'LA'));
end

% Reorganize arrays like UDG
Minv = reshape(Minv, [nd*nd np ne]);
M    = reshape(M,    [nd*nd np ne]);

Minv = permute(Minv, [2 1 3]);
M    = permute(M,    [2 1 3]);
Smin = reshape(Smin,[np 1 ne]);
Smax = reshape(Smax,[np 1 ne]);

% % Smooth the field hloc by using a CG projection
% if smoothing==1
%     [~,cgelcon,rowent2elem,colent2elem,~] = mkcgent2dgent(dgnodes,1e-6);
%     Smin = dg2cg2(Smin, cgelcon, colent2elem, rowent2elem);
% end

% % Smooth the fields M and Minv using a CG progection (2D so far)
% M(:,1,:) = dg2cg2(M(:,1,:), cgelcon, colent2elem, rowent2elem);
% M(:,2,:) = dg2cg2(M(:,2,:), cgelcon, colent2elem, rowent2elem);
% M(:,4,:) = dg2cg2(M(:,4,:), cgelcon, colent2elem, rowent2elem);
% Minv(:,1,:) = dg2cg2(Minv(:,1,:), cgelcon, colent2elem, rowent2elem);
% Minv(:,2,:) = dg2cg2(Minv(:,2,:), cgelcon, colent2elem, rowent2elem);
% Minv(:,4,:) = dg2cg2(Minv(:,4,:), cgelcon, colent2elem, rowent2elem);
% M(:,3,:)    = M(:,2,:);
% Minv(:,3,:) = Minv(:,2,:);

end

function [Jdgn, dgn] = voljac(shapgeomvt,dgnodes,smoothing)
% VOLJAC computes physical nodes, Jacobian inverse at Gauss points
%
% SYNTAX : [Jg, pg] = voljac(shapgeomvt,dgnodes)
%
% INPUTS:
%    SHAPGEOMVT : Shape functions and derivatives at DG nodes
%    DGNODES    : Geometry DG nodes
% OUTPUTS:
%    JDGN : Jacobian mapping at DG nodes
%    DGN  : Physical nodes at DG nodes (should be the same as DGNODES)
%

ngv = size(shapgeomvt,1);
npv = size(shapgeomvt,2);
nd  = size(dgnodes,2);
ne  = size(dgnodes,3);
dshapvt = reshape(permute(shapgeomvt(:,:,2:nd+1),[1 3 2]),[ngv*nd npv]);
shapgeomvt = shapgeomvt(:,:,1);
% compute dgnodes positions
pn  = reshape(dgnodes,[npv nd*ne]);
dgn = shapgeomvt*pn;
dgn = reshape(dgn,[ngv nd ne]);
% compute the Jacobian matrix of element mapping at dgnodes
Jdgn = dshapvt*pn;

for k=1:smoothing
% smoothing the jacobian mapping
    Jdgn = reshape(Jdgn,[ngv nd*nd ne]);
    [~,cgelcon,rowent2elem,colent2elem,~] = mkcgent2dgent(dgnodes,1e-6);
    for j = 1:nd*nd
        Jdgn(:,j,:) = dg2cg2(Jdgn(:,j,:), cgelcon, colent2elem, rowent2elem);
    end
end

Jdgn = permute(reshape(Jdgn,[ngv nd nd ne]),[2 3 1 4]);
Jdgn = reshape(Jdgn,[nd,nd,ngv*ne]);

end

