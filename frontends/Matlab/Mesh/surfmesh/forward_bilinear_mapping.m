function [x, y] = forward_bilinear_mapping(xi, eta, c)
% Inputs:
% xi, eta - Reference coordinates
% c - Coefficients matrix of size (4x2)
% Outputs:
% x, y - Physical coordinates


xy = [ones(length(xi(:)),1) xi(:) eta(:) xi(:).*eta(:)]*c;
x = xy(:,1);
y = xy(:,2);