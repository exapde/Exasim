
function mesh = getMeshNearBoundary(mesh,numLayers,bouIndex,elemtype,bndexpr)
%
%   DESCRIPTION:
%       Get submesh formed by "numLayers" layers of elements next ot the
%       boundary "bouIndex".
%
%   INPUTS:
%       MESH:   Original mesh
%       NUMLAYERS:  Number of layers to collect.
%       BOUINDEX: Index of the boundary whose adjacent layers of elements
%           are to be retained.
%
%   OUTPUTS:
%       MESH: Resulting mesh formed by numLayers of elements next to
%           bouIndex boundary.
%
% Written by: P. Fernandez

% Find vertices on boundary:
pf = mesh.f(mesh.f(:,4)==-abs(bouIndex),1:2);
pf = unique(pf(:));

% Find elements possessing at least 1 of these points:
ef = find(sum(ismember(mesh.t,pf),2)>0);

% Expand to multiple layers:
for i=1:numLayers-1
    % Find points encompassing one level further out:
    pf = mesh.t(ef,:);
    pf = unique(pf(:));
    
    % Find elements with those points:
    ef = find(sum(ismember(mesh.t,pf),2)>0);
end

% Update mesh:
mesh.t = mesh.t(ef,:);
mesh.dgnodes = mesh.dgnodes(:,:,ef);
[mesh.f,mesh.t2f] = mkt2f(mesh.t,elemtype);
mesh.f = setbndnbrs(mesh.p,mesh.f,bndexpr);

end
