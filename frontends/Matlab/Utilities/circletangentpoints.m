function [T1, T2] = circletangentpoints(xC, yC, r, x, y)
% CIRCLE_TANGENT_POINTS  Tangency points from an external point to a circle.
%
% Inputs:
%   (xC,yC) : circle center
%   r       : radius (positive)
%   (x,y)   : external point
%
% Outputs:
%   T1, T2  : 1x2 vectors [xT yT] for the two tangent points
%
% Notes:
%   - Requires the point to be strictly outside the circle: d > r.
%   - If d == r, there is one tangent point (the point itself on the circle).
%   - If d < r, no real tangents.

    % Vector from center to point
    vx = x - xC;
    vy = y - yC;
    d2 = vx*vx + vy*vy;
    d  = sqrt(d2);

    if r <= 0
        error('Radius r must be positive.');
    end

    if d < r
        error('Point is inside the circle: no real tangents (d < r).');
    elseif abs(d - r) < 1e-14
        % Point is on the circle: one tangent point (degenerate)
        T1 = [x, y];
        T2 = [x, y];
        return;
    end

    % Unit vector from center to point
    ux = vx / d;
    uy = vy / d;

    % Along- and perpendicular-components for tangent points
    % Tangent points in basis {u, u_perp}:
    %   T = C + r*( (r/d) u ± sqrt(1-(r/d)^2) u_perp )
    a = r / d;
    b = sqrt(1 - a*a);

    % Perpendicular unit vector (rotate u by +90°)
    px = -uy;
    py =  ux;

    % Tangent points
    tx1 = xC + r*(a*ux + b*px);
    ty1 = yC + r*(a*uy + b*py);

    tx2 = xC + r*(a*ux - b*px);
    ty2 = yC + r*(a*uy - b*py);

    T1 = [tx1, ty1];
    T2 = [tx2, ty2];
end