function [p2,t2] = removeHangingNodes(p,t,pb)
% Function [p2,t2] = removeHangingNodes(p,t,pb)
% Removes hanging nodes and associated triangles, then remeshes the 
% resulting void using Delaunay (making sure new triangles all have
% midpoints within the void created).
%
% If ordered end-repeated boundary points pb are provided, also clears
% hanging nodes from the boundaries.
%
% NOTE: this code assumes only convex voids are created (and re-meshed).
% In weird pathological cases, non-convex voids could be created, and
% this use of the Delaunay algorithm would then introduce undesired
% triangles that overlap existing ones (i.e. triangles would be added
% that have midpoints outside the void). This code can be extended
% to explicitly ensure that Delaunay's new triangles are all inside
% the void, if this proves to be a problem.
%
%--------------------------------------------------------------------------
% REVISION HISTORY:
% When     Who               What
% 08Sep13  Hemant Chaurasia  Created
%--------------------------------------------------------------------------

if nargin==2
    ih = findHangingNodes(p,t);
elseif nargin==3
    ih = findHangingNodes(p,t,pb);
else
    error('Incorrect number of input arguments.');
end
[p2,t2] = killNodes(p,t,ih);