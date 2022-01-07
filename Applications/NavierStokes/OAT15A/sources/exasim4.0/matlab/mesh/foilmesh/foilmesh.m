function mesh = foilmesh(porder,elemtype,bndexpr,xfl,yfl,xfu,yfu,R1,R2,L,nl,nu,nr1,nr2,nw,nt,a,b1,b2,c,curvedElements)
% porder    : polynomial degree
% elemtype  : 0 tri and 1 quad
% bndexpr   : expression of the boundary
% (xfl,yfl) : coordinates of the lower surface
% (xfu,yfu) : coordinates of the upper surface
% R1        : radius of the first subdomain
% R2        : radius of the whole domain
% L         : distance from the trailing edge to the far field wake
% nl        : number of edges along the lower surface 
% nu        : number of edges along the upper surface 
% nr1       : number of edges along the radial direction [0, R1]
% nr2       : number of edges along the radial direction [R1, R2]
% nw        : number of edges along the wake [0, L]
% nt        : number of edges along the open trailing edge
% a         : scaling factor along the airfoil surface     
% b1        : scaling factor along the radial direction [0, R1]
% b2        : scaling factor along the radial direction [R1, R2]
% c         : scaling factor along the wake [0, L]
% curvedElements: 0 for linear elements, 1 for curved elements
% Example   : bndexpr = {'sqrt((p(:,1)-.5).^2+p(:,2).^2)<2','true'}; 
%             mesh = foilmesh(porder,elemtype,bndexpr,xfl,yfl,xfu,yfu,0.45,10,10,100,100,50,15,40,5,1.5,4.0,4.0,4.0,1);
if nargin < 16; curvedElements = 21; end

[X12,Y12,X3,Y3] = foilcart(xfl,yfl,xfu,yfu,porder,nl,nu,nr1,nr2,nw,nt,a,b1,b2,c,R1,R2,L);

[p2,t2,dgn2]=cart2dg(elemtype,porder,X12,Y12);
[p2,t2]=fixmesh(p2,t2);

%figure(2); clf; simpplot(p2,t2); 

[p3,t3,dgn3]=cart2dg(elemtype,porder,X3,Y3);
[p3,t3]=fixmesh(p3,t3);

%figure(3); clf; simpplot(p3,t3);

[p,t] = connectmesh(p2,t2,p3,t3,1e-8);
[p,t]=fixmesh(p,t);

nodetype = 0;
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);

if curvedElements == 1
    mesh.dgnodes = cat(3,dgn2,dgn3);
end
