

function [p,t] = pt4sidesDriven(bottom,top,lateralPattern,elemtype)
%
%   DESCRIPTION:
%       Create (p,t) for 4-sides, structured, linear mesh made of triangles
%       or quadrilateral elements. The bottom and top sides can be curved.
%       The other two sides are straight lines and the node distribution
%       along them is specified through lateralPattern
%       In parametric coordinates, iso-lateral are curves whereas
%       iso-bottom/iso-top are curves.
%
%   INPUTS:
%       BOTTOM: [numNodesInBottomSide,2] (x,y) coordinates of nodes in
%               bottom side.
%       TOP:    [numNodesInTopSide,2] (x,y) coordinates of nodes in
%               top side. numNodesInTopSide = numNodesInBottomSide.
%       LATERALPATTERN: [numNodesInLateralSides] Relative position of nodes
%               in lateral sides with respect to the length of the sides
%               themselves (from 0 to 1).
%       ELEMTYPE: 0 if simplex elements (triangles). Default value if no
%                   input provided.
%                 1 if tensor product elements (quadrilaterals).
%
%   OUTPUTS:
%       P: [numNodesInMesh,2]:  (x,y) coordinates of nodes in resulting
%           mesh.
%       T: [numElementsInMesh,numVerticesPerElem]
%
%   SKETCH:
%
%                      Top
%         ___        ________________\
%       /|\  \______/               /|\
%        |                           |
%        |                           |
%  Lat 1 |                    _      | Lat 2
%        |                   / \     |
%        |     ____         /   \    |
%        |    /    \_______/     \   |
%        |___/                    \_\|
%                    Bottom         /
%
%
%   Written by: P. Fernandez
%

if nargin < 4; elemtype = 0; end

if size(bottom,1) ~= size(top,1); error('Curved sides do not have the same number of nodes.'); end

xDiv = size(bottom,1)-1;
yDiv = length(lateralPattern)-1;

% Normalize lateralPattern:
lateralPattern = (lateralPattern - lateralPattern(1))/(lateralPattern(end) - lateralPattern(1));
lateralNodes = length(lateralPattern);

xBottom = bottom(:,1);  yBottom = bottom(:,2);
xTop = top(:,1);        yTop = top(:,2);

% Compute (x,y) coordinates of vertices in mesh:
x2dArray = kron(xBottom',ones(lateralNodes,1)) + kron((xTop-xBottom)',lateralPattern(:));
y2dArray = kron(yBottom',ones(lateralNodes,1)) + kron((yTop-yBottom)',lateralPattern(:));

% Order mesh in direction of least elements for smaller bandwidth in
% FEM matrix:
if xDiv < yDiv
    minDiv = xDiv+1;
    x2dArray = permute(x2dArray,[2,1]);
    y2dArray = permute(y2dArray,[2,1]);
else
    minDiv = yDiv+1;
    % Flip to obtain a counterclockwise numbering of the vertices in t:
    x2dArray = flip(x2dArray,1);
    y2dArray = flip(y2dArray,1);
end

x1dArray = reshape(x2dArray,[],1);
y1dArray = reshape(y2dArray,[],1);

% Compute p:
p = [x1dArray,y1dArray];

% Compute t (counterclockwise ordering of vertices in element):
numElements = xDiv*yDiv;
elementIndex = (1:numElements)';
iIndex = ceil((elementIndex-0.5)/(minDiv-1));
jIndex = elementIndex - (iIndex-1)*(minDiv-1);

if elemtype == 1
    t = [(iIndex-1)*minDiv+jIndex , (iIndex-1)*minDiv+jIndex+1 , iIndex*minDiv+jIndex+1 , iIndex*minDiv+jIndex];

elseif elemtype == 0
    t1 = [(iIndex-1)*minDiv+jIndex , iIndex*minDiv+jIndex+1, iIndex*minDiv+jIndex];
    t2 = [(iIndex-1)*minDiv+jIndex , (iIndex-1)*minDiv+jIndex+1 , iIndex*minDiv+jIndex+1];
    
    t = kron(t1,[1;0]) + kron(t2,[0;1]);
end

% It may be necessary to flip t to obtain a counterclockwise numbering:
if (xBottom(2)-xBottom(1))*(yTop(1)-yBottom(1)) - (xTop(1)-xBottom(1))*(yBottom(2)-yBottom(1)) < 0
    t = flip(t,2);
end

end
