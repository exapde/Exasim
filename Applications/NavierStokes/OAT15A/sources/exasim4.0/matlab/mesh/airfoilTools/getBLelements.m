function mesh = getBLelements(mesh,NLAYERS,XLIM,bndexpr,PLOT)
% Function to find a layer of elements adjacent to an airfoil boundary
% (assumed to be boundary #1).
%
%   mesh = getBLelements(mesh,NLAYERS,XLIM,fd,PLOT)
%
%     mesh    = C mesh
%     NLAYERS = number of extrusion layers to preserve
%     XLIM    = upper bound on x coordinate for preserved elements
%               OR, could also be a pair of points defining a line.
%     bndexpr = boundary expressions
%     PLOT    = plotting flag
%
%--------------------------------------------------------------------------
% REVISION HISTORY:
% When     Who               What
% 03May13  Hemant Chaurasia  Created
% 18Jun13  Hemant Chaurasia  Begin update to take XLIM as a pair of points.
%--------------------------------------------------------------------------
if nargin<5; PLOT=0; end
% clear;
% mesh=mkmesh_naca0006(4,0,2.1,0);
% NLAYERS = 4;
% XLIM = 1.125;

% Find points on airfoil boundary
pf = mesh.f(mesh.f(:,4)==-1,1:2);
pf = unique(pf(:));
if PLOT; plot(mesh.p(pf,1),mesh.p(pf,2),'o'); axis equal; end

% Find elements possessing at least 1 of these points
ef = find(sum(ismember(mesh.t,pf),2)>0);

% Expand to multiple layers
for i=1:NLAYERS-1
    pf = mesh.t(ef,:); pf = unique(pf(:));   % Find points encompassing one level further out
    ef = find(sum(ismember(mesh.t,pf),2)>0); % Find elements with those points
end

% Restrict to x<XLIM
if size(XLIM,1)==1 && size(XLIM,2)==1
    ef = ef(max(mesh.dgnodes(:,1,ef))<XLIM);
elseif size(XLIM,1)==2 && size(XLIM,2)==2
    error('Not yet implemented.');
else
    error('Did not recognize XLIM.');
end

% Construct submesh:
mesh.t = mesh.t(ef,:);
mesh.dgnodes = mesh.dgnodes(:,:,ef);
elemtype = 0;
[mesh.f,mesh.t2f] = mkt2f(mesh.t,elemtype);
mesh.f = setbndnbrs(mesh.p,mesh.f,bndexpr);

% % Diagnostic plots
% plotbndfaces(mesh,1);
% figure; hold on; plot(mesh.p(pf,1),mesh.p(pf,2),'bo');
% plot(x,y,'rx'); hold off;

% i1 = findbndpoints(mesh,1);
% i2 = findbndpoints(mesh,2);
% figure; hold on;
% plot(mesh.p(i1,1),mesh.p(i1,2),'bo');
% plot(mesh.p(i2,1),mesh.p(i2,2),'rx');
% hold off;

if PLOT
    meshplot(mesh,1);
    disp(['BL contains ' num2str(length(ef)) ' elements.']);
end