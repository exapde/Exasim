
function [x,y] = cmeshTransition( nxw, nflr, nflf, nfuf, nfur, nr, sps, spr, transitionLocation, numRadialLayers, numAxialElementsWake)
%CMESHPARAM  Creates mesh in parametric space for airfoil c-type grids
%
% x expands from -2 to 2. Airfoil expands from -1 to 1.
% y expands from 0 to 1. This makes it easie to make C-mesh "socks" for 
%       embedding in distmesh farfield meshes (compared to yMax = 0.5)
%        __________________________________ ______
%       |      |      |      |      |      |      |
%       |      |      |      |      |      |      |
%    nr |      |      |      |      |      |      |
%       |      |      |      |      |      |      |
%       |______|______|______|______|______|______|
%          nxw   nflr   nflf   nfuf   nfur
%
%   nxw  : number of subdivison in the wake
%   nflr : number of subdivision in the lower foil (after transition)
%   nflf : number of subdivision in the lower foil (before transition)
%   nfuf : number of subdivisions in the upper foil (before transition)
%   nfur : number of subdivisions in the upper foil (after transition)
%   nr   : number of subdivisions in the radial direction
%   sps(id) : streamwise size control
%     sps(1)  - ratio between the first and last elements in the wake 
%     sps(2)  - ratio between transition point and trailing edge element (lower)
%     sps(3)  - ratio between transition point and leading edge element (lower)
%     sps(4)  - ratio between transition point and leading edge element (upper)
%     sps(5)  - ratio between transition point and trailing edge element (upper)
%     sps(6)  - ratio between the first and last elements in the wake (lower far field)
%     sps(7)  - ratio between transition point and trailing edge element (lower far field)
%     sps(8)  - ratio between transition point and leading edge element (lower far field)
%     sps(9)  - ratio between transition point and leading edge element (upper far field)
%     sps(10) - ratio between transition point and trailing edge element (upper far field)
%     sps(11) - ratio between the first and last elements in the wake (upper far field)
%   sps(id) : radial size control
%     spr(1) - ratio between far-field and wake element (lower) 
%     spr(2) - ratio between far-field and trailing edge (lower)
%     spr(3) - ratio between far-field and transition point (lower)
%     spr(4) - ratio between far-field and leading edge
%     spr(5) - ratio between far-field and transition point (upper)
%     spr(6) - ratio between far-field and trailing edge (upper)
%     spr(7) - ratio between far-field and wake element (upper)
%   transitionLocation: Transition location (x_trans / c)
%     transitionLocation(1): Upper surface
%     transitionLocation(2): Lower surface

if nargin < 9; transitionLocation = [0.5,0.5]; end
if nargin < 10; numRadialLayers = nr; end
if nargin < 11; numAxialElementsWake = nxw; end

if sum(transitionLocation < 0) || sum(transitionLocation > 1)
    error('Transition location must be in range [0,1].')
end

YED = 1.0;
yMax = 1.0;
wakeExtension = 1.0;

X1 = block([ nxw,nr],[1/sps(1),  1/sps(6), spr(1), spr(2)], [-1.0-wakeExtension, 0.0; -1.0, 0.0; -1.0-wakeExtension, YED; -1.0, yMax]);
X2 = block([nflr,nr],[  sps(2),    sps(7), spr(2), spr(3)], [-1.0, 0.0; -transitionLocation(2), 0.0; -1.0, yMax; -transitionLocation(2), yMax]);
X3 = block([nflf,nr],[1/sps(3),  1/sps(8), spr(3), spr(4)], [-transitionLocation(2), 0.0;  0.0, 0.0; -transitionLocation(2), yMax;  0.0, yMax]);
X4 = block([nfuf,nr],[  sps(4),    sps(9), spr(4), spr(5)], [ 0.0, 0.0;  transitionLocation(1), 0.0;  0.0, yMax;  transitionLocation(1), yMax]);
X5 = block([nfur,nr],[1/sps(5), 1/sps(10), spr(5), spr(6)], [ transitionLocation(1), 0.0;  1.0, 0.0;  transitionLocation(1), yMax;  1.0, yMax]);
X6 = block([ nxw,nr],[  sps(1),   sps(11), spr(6), spr(7)], [ 1.0, 0.0;  1.0+wakeExtension, 0.0;  1.0, yMax;  1.0+wakeExtension, YED]);

X1 = X1(:,end-numAxialElementsWake+1:end,1:numRadialLayers);
X2 = X2(:,:,1:numRadialLayers);
X3 = X3(:,:,1:numRadialLayers);
X4 = X4(:,:,1:numRadialLayers);
X5 = X5(:,:,1:numRadialLayers);
X6 = X6(:,1:numAxialElementsWake,1:numRadialLayers);

x = [squeeze(X1(1,:,:))', squeeze(X2(1,2:end,:))', squeeze(X3(1,2:end,:))', ...
                          squeeze(X4(1,2:end,:))', squeeze(X5(1,2:end,:))', squeeze(X6(1,2:end,:))'];
y = [squeeze(X1(2,:,:))', squeeze(X2(2,2:end,:))', squeeze(X3(2,2:end,:))', ...
                          squeeze(X4(2,2:end,:))', squeeze(X5(2,2:end,:))', squeeze(X6(2,2:end,:))'];

end
