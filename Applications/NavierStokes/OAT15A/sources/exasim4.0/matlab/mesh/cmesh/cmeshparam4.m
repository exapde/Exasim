function [x,y] = cmeshparam( nxw, nfl, nfu, nr, sps, spr)
%CMESHPARAM  Creates mesh in parametric space for airfoil c-type grids
%       _______________________________
%       |       |       |       |       |
%       |       |       |       |       |
%    nr |       |       |       |       |
%       |       |       |       |       |
%       |_______|_______|_______|_______|
%          nxw     nfl     nfu
%
%   nxw : number of subdivison in the wake
%   nfl : number of subdivision in the lower foil 
%   nfu : number of subdivisions in the upper foil
%   nr  : number of subdivisions in the radial direction
%   sps(id) : streamwise size control
%     sps(1) - ratio between the first and last elements in the wake 
%     sps(2) - ratio between leading edge and trailing edge element (lower)
%     sps(3) - ratio between leading edge and trailing edge element (upper)
%     sps(4) - ratio between the first and last elements in the wake (lower far field)
%     sps(5) - ratio between leading edge and trailing edge element (lower far field)
%     sps(6) - ratio between leading edge and trailing edge element (upper far field)
%     sps(7) - ratio between the first and last elements in the wake (upper far field)
%   sps(id) : radial size control
%     spr(1) - ratio between far-field and wake element (lower) 
%     spr(2) - ratio between far-field and trailing edge (lower)
%     spr(3) - ratio between far-field and leading edge
%     spr(4) - ratio between far-field and trailing edge (upper)
%     spr(5) - ratio between far-field and wake element (upper)

X1 = block([nxw,nr],[1/sps(1), 1/sps(4), spr(1), spr(2)], [-2,0;-1,0;-2,0.5;-1,1]);
X2 = block([nfl,nr],[  sps(2),   sps(5), spr(2), spr(3)], [-1,0; 0,0;-1,1; 0,1]);
X3 = block([nfu,nr],[1/sps(3), 1/sps(6), spr(3), spr(4)], [ 0,0; 1,0; 0,1; 1,1]);
X4 = block([nxw,nr],[  sps(1),   sps(7), spr(4), spr(5)], [ 1,0; 2,0; 1,1; 2,0.5]);

% X1 = block([nxw,nr],[1/sps(1), 1/sps(4), spr(1), spr(2)], [-2,0;-1,0;-2,1;-1,1]);
% X2 = block([nfl,nr],[  sps(2),   sps(5), spr(2), spr(3)], [-1,0; 0,0;-1,1; 0,1]);
% X3 = block([nfu,nr],[1/sps(3), 1/sps(6), spr(3), spr(4)], [ 0,0; 1,0; 0,1; 1,1]);
% X4 = block([nxw,nr],[  sps(1),   sps(7), spr(4), spr(5)], [ 1,0; 2,0; 1,1; 2,1]);

x = [squeeze(X1(1,:,:))', squeeze(X2(1,2:end,:))', squeeze(X3(1,2:end,:))', squeeze(X4(1,2:end,:))'];
y = [squeeze(X1(2,:,:))', squeeze(X2(2,2:end,:))', squeeze(X3(2,2:end,:))', squeeze(X4(2,2:end,:))'];
