function [x,y] = cmeshparam( nxw, nflr, nflf, nfuf, nfur, nr, sps, spr)
%CMESHPARAM  Creates mesh in parametric space for airfoil c-type grids
%        __________________________________ ______
%       |      |      |      |      |      |      |
%       |      |      |      |      |      |      |
%    nr |      |      |      |      |      |      |
%       |      |      |      |      |      |      |
%       |______|______|______|______|______|______|
%          nxw   nflr   nflf   nfuf   nfur
%
%   nxw  : number of subdivison in the wake
%   nflr : number of subdivision in the lower foil (rear)
%   nflf : number of subdivision in the lower foil (front)
%   nfuf : number of subdivisions in the upper foil (front)
%   nfur : number of subdivisions in the upper foil (rear)
%   nr   : number of subdivisions in the radial direction
%   sps(id) : streamwise size control
%     sps(1)  - ratio between the first and last elements in the wake 
%     sps(2)  - ratio between mid-chord and trailing edge element (lower)
%     sps(3)  - ratio between mid-chord and leading edge element (lower)
%     sps(4)  - ratio between mid-chord and leading edge element (upper)
%     sps(5)  - ratio between mid-chord and trailing edge element (upper)
%     sps(6)  - ratio between the first and last elements in the wake (lower far field)
%     sps(7)  - ratio between mid-chord and trailing edge element (lower far field)
%     sps(8)  - ratio between mid-chord and leading edge element (lower far field)
%     sps(9)  - ratio between mid-chord and leading edge element (upper far field)
%     sps(10) - ratio between mid-chord and trailing edge element (upper far field)
%     sps(11) - ratio between the first and last elements in the wake (upper far field)
%   sps(id) : radial size control
%     spr(1) - ratio between far-field and wake element (lower) 
%     spr(2) - ratio between far-field and trailing edge (lower)
%     spr(3) - ratio between far-field and mid-chord (lower)
%     spr(4) - ratio between far-field and leading edge
%     spr(5) - ratio between far-field and mid-chord (upper)
%     spr(6) - ratio between far-field and trailing edge (upper)
%     spr(7) - ratio between far-field and wake element (upper)

X1 = block([nxw,nr],[1/sps(1),  1/sps(6), spr(1), spr(2)], [-2.0, 0.0; -1.0, 0.0; -2.0, 0.5; -1.0, 1.0]);
X2 = block([nxw,nr],[  sps(2),    sps(7), spr(2), spr(3)], [-1.0, 0.0; -0.5, 0.0; -1.0, 1.0; -0.5, 1.0]);
X3 = block([nxw,nr],[1/sps(3),  1/sps(8), spr(3), spr(4)], [-0.5, 0.0;  0.0, 0.0; -0.5, 1.0;  0.0, 1.0]);
X4 = block([nfl,nr],[  sps(4),    sps(9), spr(4), spr(5)], [ 0.0, 0.0;  0.5, 0.0;  0.0, 1.0;  0.5, 1.0]);
X5 = block([nfl,nr],[1/sps(5), 1/sps(10), spr(5), spr(6)], [ 0.5, 0.0;  1.0, 0.0;  0.5, 1.0;  1.0, 0.5]);
X6 = block([nfu,nr],[  sps(1),   sps(11), spr(6), spr(7)], [ 1.0, 0.0;  2.0, 0,0;  1.0, 0.5;  2.0, 0,0]);


% X1 = block([nxw,nr],[1/sps(1), 1/sps(4), spr(1), spr(2)], [-2,0;-1,0;-2,1;-1,1]);
% X2 = block([nfl,nr],[  sps(2),   sps(5), spr(2), spr(3)], [-1,0; 0,0;-1,1; 0,1]);
% X3 = block([nfu,nr],[1/sps(3), 1/sps(6), spr(3), spr(4)], [ 0,0; 1,0; 0,1; 1,1]);
% X4 = block([nxw,nr],[  sps(1),   sps(7), spr(4), spr(5)], [ 1,0; 2,0; 1,1; 2,1]);

x = [squeeze(X1(1,:,:))', squeeze(X2(1,2:end,:))', squeeze(X3(1,2:end,:))', ...
                          squeeze(X4(1,2:end,:))', squeeze(X5(1,2:end,:))', squeeze(X6(1,2:end,:))'];
y = [squeeze(X1(2,:,:))', squeeze(X2(2,2:end,:))', squeeze(X3(2,2:end,:))', ...
                          squeeze(X4(2,2:end,:))', squeeze(X4(2,2:end,:))', squeeze(X4(2,2:end,:))'];
