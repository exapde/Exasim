function [x,y] = hmeshparam %( nxw, nfl, nfu, nr, sps, spr)
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

X1 = block([10,10],[1/10,1/3,3,15], [-1,0; 0,0;-1,1; 0,1]);
X2 = block([20,10],[-8,1,15,15], [ 0,0; 1,0; 0,1; 1,1]);
X3 = block([10,10],[10,3,15,3], [ 1,0; 2,0; 1,1; 2,1]);

x = [squeeze(X1(1,:,:))', squeeze(X2(1,2:end,:))', squeeze(X3(1,2:end,:))'];
y = [squeeze(X1(2,:,:))', squeeze(X2(2,2:end,:))', squeeze(X3(2,2:end,:))'];
