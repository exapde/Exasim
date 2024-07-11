function mesh = mkmesh_gap(m,n,o,porder,a,b,c,delta,elemtype,nodetype)
%MKMESH_GAP Creates 2 3D meshes separated by a gap
% The mesh created is composed of 2 parts. One bottom, and one top part
% separated by a gap delta at height z \approx c/2. The interface has a
% sinusoidal shape.
%
% SYNTAX :   mesh = mkmesh_gap(m,n,o,porder,a,b,c,delta,elemtype,nodetype)
%
% INPUTS:
%      m:         Number of points in the x direction 
%      n:         Number of points in the y direction
%      o:         Number of points in the z direction
%      porder:    Polynomial Order of Approximation (default=1)
%      a:         Length in the x direction
%      b:         Length in the y direction
%      c:         Length in the z direction
%      delta:     length of the gap
%      elemtype:  Flag determining element type
%                 Flag = 0 tri/tet elements (default)
%                 Flag = 1 quad/hex elements
%      nodetype:  Flag determining node distribution 
%                 Flag = 0 uniform distribution (default)
%                 Flag = 1 nonuniform distribution
%
% OUTPUTS:
%      mesh:      Mesh structure
%
%   See also: CUBEMESH, MKMESH, MKMESH_CUBE
%
% Author(s): Sebastien TERRANA
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email address: XXXXXXX@mit.edu 
% Website: http://aeroastro.mit.edu/
% March 2018

if nargin<1, m=2; end
if nargin<2, n=m; end
if nargin<3, o=n; end
if nargin<4, porder = 1; end


if m < 2 || n < 2 || o < 2,
    error('At least m=2, n=2, o=2 needed.');
end

% Boundary Expression for mkmesh
bndexpr = {'all(p(:,1)<min(p0(:,1))+1e-6)','all(p(:,1)>max(p0(:,1))-1e-6)', ...
           'all(p(:,2)<min(p0(:,2))+1e-6)','all(p(:,2)>max(p0(:,2))-1e-6)', ...
           'all(p(:,3)<min(p0(:,3))+1e-6)','all(p(:,3)>max(p0(:,3))-1e-6)', ...
           'all(abs(p(:,3)+1)<1e-6)','all(abs(p(:,3))<1e-6)','true'};  

% Nb of pts in z direction for bottom and top meshes
o  = o + 1;
o1 = floor(o/2); o2 = o-o1;

% Bottom part of the mesh
[p,t] = cubemesh(m,n,o1,elemtype);
np1 = size(p,1); ne1 = size(t,1);
% Top part of the mesh
[p2,t2] = cubemesh(m,n,o2,elemtype);

% Gathering the two meshes (meshes are separated in case of node collapsed)
p(:,3) = p(:,3)-2;
p  = [p;p2];
t2 = t2+np1;
t  = [t;t2];

% Creating the DG mesh
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);

% Get first mesh back in its initial position
mesh.p(1:ne1,3) = mesh.p(1:ne1,3) + 2;
mesh.dgnodes(:,3,1:ne1) = mesh.dgnodes(:,3,1:ne1) + 2;

% Sinus Deformation of meshes
amp = 0.05;
omega = 4*pi;
omega = 0; % Flat interface
Dgap = delta/c;

% up deformation of the bottom mesh
blk = p(1:ne1,:);
blk(:,3) = blk(:,3) + blk(:,3) .* (amp*sin(omega*blk(:,1)).*sin(omega*blk(:,2)) - Dgap);
p(1:ne1,3) = blk(:,3);
blk = mesh.dgnodes(:,1:3,1:ne1);
blk(:,3,:) =  blk(:,3,:) +  blk(:,3,:) .* (amp*sin(omega*blk(:,1,:)).*sin(omega*blk(:,2,:)) - Dgap);
mesh.dgnodes(:,3,1:ne1) = blk(:,3,:);

% down deformation of the up mesh
blk = p(ne1+1:end,:);
blk(:,3) = blk(:,3) + (1-blk(:,3)) .* (amp*sin(omega*blk(:,1)).*sin(omega*blk(:,2)) + Dgap);
p(ne1+1:end,3) = blk(:,3);
blk = mesh.dgnodes(:,1:3,ne1+1:end);
blk(:,3,:) =  blk(:,3,:) + (1-blk(:,3,:)) .* (amp*sin(omega*blk(:,1,:)).*sin(omega*blk(:,2,:)) + Dgap);
mesh.dgnodes(:,3,ne1+1:end) = blk(:,3,:);

% deformation of vertices
p(:,1) = a*p(:,1);
p(:,2) = b*p(:,2);
p(:,3) = c/2*p(:,3);
p(ne1+1:end,3) = delta + c/2 + p(ne1+1:end,3);
% deformation of HDG nodes
mesh.dgnodes(:,1,:) = a*mesh.dgnodes(:,1,:);
mesh.dgnodes(:,2,:) = b*mesh.dgnodes(:,2,:);
mesh.dgnodes(:,3,:) = c/2*mesh.dgnodes(:,3,:);
mesh.dgnodes(:,3,ne1+1:end) = c/2 + mesh.dgnodes(:,3,ne1+1:end);

