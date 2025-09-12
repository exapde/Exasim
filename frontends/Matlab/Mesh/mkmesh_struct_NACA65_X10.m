
function mesh = mkmesh_struct_NACA65_X10(X,porder,C,a,plotFlag)

% TO BE CLEANED.

% Default parameters provide nice mesh for p=4, C=0.7

if nargin<5; plotFlag = 0; end
if nargin<4; a = 1.0; end
if nargin<3; C = 0.7; end           % Mesh scaling parameter
if nargin<2; porder = 4; end

hybrid = 'hdg';
BP = 1;
chord = 1.0;                                                        % chord length scaling
alpha = (pi/180)*3.5;               % Angle of attack.
theta = (pi/180)*45 - alpha; % atan(2.5/5.0);                                                % anticlockwise rotation angle applied to airfoil
xTranslation = 0;                                                          % x-translation applied to scaled & rotated sock mesh
yTranslation = 0;                                                          % y-translation applied to scaled & rotated sock mesh
elemtype = 1;
nodetype = 0;
plotCmeshFlag = 0;

xLeftDomain = -1;                                                          % left boundary x-coord
xStartDomain2 = -0.7;
xRightDomain = 3;                                                          % right boundary x-coord
yLeftBottom = xLeftDomain*tan(theta) - BP/2;                                                       % y-coord of bottom-left corner of blade passage mesh
yLeftTop = yLeftBottom + BP;
yRightBottom = xRightDomain*tan(theta) - BP/2;   
yRightTop = yRightBottom + BP;

nLeft = round(C*18);                                                % number of intervals on left boundary
nTop1 = round(C*30);  % number of intervals on top boundary, up to X1
nTop2 = round(C*30);  % number of intervals on top boundary, between X1 and X2
nTop3 = round(C*48);  % number of intervals on top boundary, past X2
nRight = round(C*30);                                                % number of intervals on right boundary
nBottom3 = round(C*48);  % number of intervals on top boundary, past X2
nBottom2 = round(C*30);  % number of intervals on top boundary, between X1 and X2
nBottom1 = round(C*30);  % number of intervals on top boundary, up to X1

transitionLocation = [0.5,0.5];

% C-mesh sock parameters:
cmeshparams.n1  = round(C*60)*porder+1;                                   % Number of subdivison in the wake.
cmeshparams.n2  = round(C*38)*porder+1;                                   % Number of subdivision in the front part of the foil (both upper and lower surfaces).
cmeshparams.n2r = round(C*75)*porder+1;                                   % Number of subdivision in the front part of the foil (both upper and lower surfaces).
cmeshparams.n3  = round(C*50)*porder+1;                                   % Number of subdivisions in the radial direction.
TEC       = 5;                                                      % TE Compression
cmeshparams.sps = [TEC, TEC, 1, 1, TEC, TEC, 1, 1, 1, 1, TEC];            % Streamwise size control.
cmeshparams.spr = [10, 10, 10, 10, 10, 10, 10]*30;                        % Ratio of cell size in the radial direction between the far field and the closest cell to the foil.

numLayersInCmesh = round(12*C);                    % number of layers of elements in C-mesh sock

airfoilBouIndex = 1; % Since we use bndexpr = {'dist_NACA65_010(p)<0.001','true'} for the cmesh.
outerBouIndexSock = 2; % Since we use bndexpr = {'dist_NACA65_010(p)<0.001','true'} for the cmesh.

divTopBottomVsLeftDomain2 = 0.80;

growthFactorTopBottom2 = 0.99;        % Growth factor in top and bottom parts of domain 2.
radialDivDomain2 = 10;
growthFactorRadialDivDomain2 = 1.23;

axialDivDomain1 = 2;
growthFactorAxialDivDomain1 = 1.2;

growthFactorAxialDivDomains345 = 1.03;
axialDivDomains345 = 40; %80;

yScalingSide2Domain3 = 5;
yTraslationSide2Domain3 = 0.01; % tan(alpha)*(xRightDomain-chord*cos(theta));

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

writeAirfoilScalingRotationAndTranslation(chord,theta,xTranslation,yTranslation)

% Make C-mesh:
bndexpr = {['dist_NACA65_X10(p,',num2str(X),',',num2str(a),')<0.001'],'true'};
cmesh = mkcmesh_NACA65_X10(X,porder,cmeshparams,elemtype,bndexpr,a,plotCmeshFlag,transitionLocation);

% Extract boundary layer part of C-mesh
sock = getMeshNearBoundary(cmesh,numLayersInCmesh,airfoilBouIndex,elemtype,bndexpr);
sock.fcurved = ones(size(sock.fcurved));
sock.tcurved = ones(size(sock.tcurved));

% Extract outer edge points of sock
ip = findbndpoints(sock,outerBouIndexSock);
ps = sock.p(ip,:);    % sock interface points, in order, with first and last points same

[~,ipSorted] = sort(ps(2:end,1),'descend');

pEndWake = ps(ipSorted(1:2*numLayersInCmesh+1)+1,:);
pEndWake = sortrows(pEndWake,2);        % Sorted from bottom to top.

ipCurved = setdiff(2:size(ps,1),ipSorted(1:2*numLayersInCmesh+1)+1);
side1 = ps(ipCurved,:);              % Sorted from upper TE to lower TE.
side1 = [pEndWake(end,:);side1;pEndWake(1,:)];      % Include first and last node (vertex nodes)

% Scale chord length:
sock.dgnodes = sock.dgnodes*chord;
sock.p = sock.p*chord;
pEndWake = pEndWake*chord;
side1 = side1*chord;

% Rotate anticlockwise:
A_ROT = [cos(theta) -sin(theta); sin(theta) cos(theta)];
for i=1:size(sock.dgnodes,3)
    sock.dgnodes(:,:,i) = (A_ROT*sock.dgnodes(:,:,i)')';
end
sock.p = (A_ROT*sock.p')';
pEndWake = (A_ROT*pEndWake')';
side1 = (A_ROT*side1')';

% Translate in (x,y) directions:
sock.dgnodes(:,1,:) = sock.dgnodes(:,1,:) + xTranslation;
sock.dgnodes(:,2,:) = sock.dgnodes(:,2,:) + yTranslation;
sock.p(:,1) = sock.p(:,1) + xTranslation;
sock.p(:,2) = sock.p(:,2) + yTranslation;
pEndWake(:,1) = pEndWake(:,1) + xTranslation;
pEndWake(:,2) = pEndWake(:,2) + yTranslation;
side1(:,1) = side1(:,1) + xTranslation;
side1(:,2) = side1(:,2) + yTranslation;

% Make side2:
numDivSide1 = size(side1,1)-1;
numDivTopBottom2 = floor(divTopBottomVsLeftDomain2*numDivSide1/2);
numDivLeft2 = numDivSide1 - 2*numDivTopBottom2;
xEndDomain2 = 0.5*(pEndWake(1,1)+pEndWake(end,1));

relativeSize_aux = growthFactorTopBottom2.^(1:numDivTopBottom2);
relativeSize_aux = [0,cumsum(relativeSize_aux)/sum(relativeSize_aux)];

xBottom = xStartDomain2 + (xEndDomain2-xStartDomain2)*relativeSize_aux;
yBottom = -BP/2 + tan(theta)*xBottom;

xTop = flip(xBottom);
yTop = BP/2 + tan(theta)*xTop;

xLeft = xTop(end)*ones(numDivLeft2+1,1);
yLeft = linspace(yTop(end),yBottom(1),numDivLeft2+1);

% Remove duplicated nodes:
xBottomReduced = xBottom(2:end);   yBottomReduced = yBottom(2:end);
xLeftReduced = xLeft(2:end);       yLeftReduced = yLeft(2:end);

side2 = [xTop(:) yTop(:); xLeftReduced(:) yLeftReduced(:); xBottomReduced(:) yBottomReduced(:)];

% Generate p,t for Domain 2:
patternRadialDivDomain2 = growthFactorRadialDivDomain2.^(1:radialDivDomain2);
patternRadialDivDomain2 = [0,cumsum(patternRadialDivDomain2)/sum(patternRadialDivDomain2)];
[po2,to2] = pt4sidesDriven(side1,side2,patternRadialDivDomain2,elemtype);

% Generate p,t for Domain 1:
patternAxialDivDomain1 = growthFactorAxialDivDomain1.^(1:axialDivDomain1);
patternAxialDivDomain1 = [0,cumsum(patternAxialDivDomain1)/sum(patternAxialDivDomain1)];

side1_domain1 = [xLeft(:) yLeft(:)];
side2_domain1 = side1_domain1;
side2_domain1(:,1) = xLeftDomain;
side2_domain1(:,2) = side1_domain1(:,2) + tan(theta)*(side2_domain1(:,1)-side1_domain1(:,1));

[po1,to1] = pt4sidesDriven(side1_domain1,side2_domain1,patternAxialDivDomain1,elemtype);

% Generate p,t for Domain 3:
patternAxialDivDomains345 = growthFactorAxialDivDomains345.^(1:axialDivDomains345);
patternAxialDivDomains345 = [0,cumsum(patternAxialDivDomains345)/sum(patternAxialDivDomains345)];

side1_domain3 = pEndWake;
side2_domain3 = side1_domain3;
side2_domain3(:,1) = xRightDomain;

ySide2_domain3_tmp = side1_domain3(:,2) + tan(theta)*(side2_domain3(:,1)-side1_domain3(:,1));
side2_domain3(:,2) = linspace(min(ySide2_domain3_tmp),max(ySide2_domain3_tmp),size(pEndWake,1));

avgYside2 = 0.5*(side2_domain3(1,2)+side2_domain3(end,2));
side2_domain3(:,2) = yTraslationSide2Domain3 + avgYside2 + yScalingSide2Domain3*(side2_domain3(:,2)-avgYside2);

[po3,to3] = pt4sidesDriven(side1_domain3,side2_domain3,patternAxialDivDomains345,elemtype);

% Generate p,t for Domain 4:
if radialDivDomain2 >= numDivSide1;
    error('Number of divisions along airfoil must be greater than radial divisions in external mesh. The generation of (p,t) for Domains 4 and 5 is based on this hypothesis');
end

side1_domain4 = po2(1:radialDivDomain2+1,:);
side2_domain4 = side1_domain4;
side2_domain4(:,1) = xRightDomain;
yBottomSide2_domain4 = max(side2_domain3(:,2));
yTopSide2_domain4 = yRightTop;
side2_domain4(:,2) = linspace(yTopSide2_domain4,yBottomSide2_domain4,radialDivDomain2+1);

[po4,to4] = pt4sidesDriven(side1_domain4,side2_domain4,patternAxialDivDomains345,elemtype);

% Generate p,t for Domain 5:
if radialDivDomain2 >= numDivSide1;
    error('Number of divisions along airfoil must be greater than radial divisions in external mesh. The generation of (p,t) for Domains 4 and 5 is based on this hypothesis');
end

side1_domain5 = po2(end-radialDivDomain2:end,:);
side2_domain5 = side1_domain5;
side2_domain5(:,1) = xRightDomain;
yTopSide2_domain5 = min(side2_domain3(:,2));
yBottomSide2_domain5 = yRightBottom;
side2_domain5(:,2) = linspace(yBottomSide2_domain5,yTopSide2_domain5,radialDivDomain2+1);

[po5,to5] = pt4sidesDriven(side1_domain5,side2_domain5,patternAxialDivDomains345,elemtype);

% Connect meshes (important to perform this in order. sock mesh must be
% first):
[pp,tt] = connectmesh(sock.p,sock.t,po2,to2);

[pp,tt] = connectmesh(pp,tt,po1,to1);

[pp,tt] = connectmesh(pp,tt,po3,to3);

[pp,tt] = connectmesh(pp,tt,po4,to4);

[pp,tt] = connectmesh(pp,tt,po5,to5);

% Build mesh data structure from pp,tt:
bndexpr = {['dist_NACA65_X10_afterScalingAndRotation(p,',num2str(X),',',num2str(a),')<0.05'],...
           ['p(:,1)<' num2str(xLeftDomain) '+1e-3'],...
           ['p(:,1)>' num2str(xRightDomain) '-1e-3'],...
           ['dsegment(p,[' num2str(xLeftDomain) ' ' num2str(yLeftBottom) ';' num2str(xRightDomain) ' ' num2str(yRightBottom) '])<1e-3'],...
           ['dsegment(p,[' num2str(xLeftDomain) ' ' num2str(yLeftTop) ';' num2str(xRightDomain) ' ' num2str(yRightTop) '])<1e-3']};
curvedBoundary = 1;
mesh = mkmesh(pp,tt,porder,bndexpr,elemtype,nodetype,curvedBoundary);
mesh = joinBoundaries(mesh,4,5,1);
mesh.fcurved = ones(size(mesh.fcurved));  % declare all faces & tris curved
mesh.tcurved = ones(size(mesh.tcurved));
mesh.bndexpr = bndexpr;

% Stamp curved sock.dgnodes onto mesh.dgnodes:
mesh.dgnodes(:,:,1:size(sock.dgnodes,3)) = sock.dgnodes;        % This relies on the fact that the first elements in the mesh are those from the sock.


% Find outer boundary of the sock and match dgnodes to outer mesh
% (preserve the sock dgnodes and displace the outer mesh dgnodes at the interface)
sock_bf = find(sock.f(:,4)==-abs(outerBouIndexSock));                % faces of outer boundary of sock
[btri,fnum] = find(ismember(sock.t2f,sock_bf)); % triangle numbers and local face numbers in sock (same as in assembled mesh)
bf = 0*btri;
for i=1:length(btri)
    bf(i) = mesh.t2f(btri(i),fnum(i));          % face numbers in the assembled mesh
end
inner_tri = []; outer_tri = []; inner_face = []; outer_face = [];
for i=1:length(bf)                              % match faces to inner and outer triangles & local face numbers
    [temp_tri,temp_face] = find(mesh.t2f==bf(i));
    if length(temp_tri)~=2
        warning('Couldn''t find two elements for this face.');
    else
        if temp_tri(1)<=size(sock.t2f,1)
            inner_tri(end+1) = temp_tri(1); inner_face(end+1) = temp_face(1);
            outer_tri(end+1) = temp_tri(2); outer_face(end+1) = temp_face(2);
        else
            inner_tri(end+1) = temp_tri(2); inner_face(end+1) = temp_face(2);
            outer_tri(end+1) = temp_tri(1); outer_face(end+1) = temp_face(1);
        end
    end
end


% Overwrite outer dgnodes with inner ones:
for i=1:length(inner_tri)
    mesh.dgnodes(flipud(mesh.perm(:,outer_face(i))),:,outer_tri(i)) ...
    = mesh.dgnodes(mesh.perm(:,inner_face(i)),:,inner_tri(i));
end

% Reorder faces:
[mesh.f,mesh.t2f,mesh.f2f,mesh.flev] = faceordering(mesh.f, mesh.t2f); 
 
mesh.bf = reshape(mesh.f(abs(mesh.t2f'),end),[size(mesh.perm,2) size(mesh.t,1)]);

if strcmp(hybrid,'edg')
    mesh.f(:,end+1) = 1;
elseif strcmp(hybrid,'hdg')    
    mesh.f(:,end+1) = 0;
elseif strcmp(hybrid,'hedg')
    a = mesh.f(:,end);
    i = (a>0);
    mesh.f(:,end+1) = 0;
    mesh.f(i,end)   = 1;
elseif ~isempty(hybrid) 
    mesh.f = feval(hybrid,mesh);
else 
    error('Hybrid flag is not valid');
end

[mesh.elcon,mesh.nsiz] = elconnectivities(mesh);
mesh.f(:,end) = [];

mesh.dgnodes(:,3:4,:) = mesh.dgnodes(:,1:2,:);
mesh.dgnodes(:,5:6,:) = 0;

if plotFlag; meshplot(mesh,[1 1 0 0 0]); end

end
