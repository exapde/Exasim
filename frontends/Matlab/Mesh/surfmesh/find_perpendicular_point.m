function [XD, YD, position] = find_perpendicular_point(XA, YA, XB, YB, XC, YC)
    % Determine if line AB is vertical
    if XB == XA
        XD = XA;
        YD = YC; % Perpendicular line passes through the same x-coordinate

        % Determine if C is to the left or right of the vertical line
        if XC > XA
            position = 1;  %'right of line AB';
        elseif XC < XA
            position = -1; % 'left of line AB';
        else
            position = 0;  %'on line AB';
        end
        return;
    end

    % Slope of line AB
    m = (YB - YA) / (XB - XA);

    % Y-intercept of line AB
    c = YA - m * XA;

    % Solve for perpendicular point D
    m_perp = -1 / m; % Perpendicular slope
    XD = ((YC - m_perp * XC) - c) / (m - m_perp);
    YD = m * XD + c;

    % Determine position of C relative to AB
    % Substitute XC into the line equation to get the y-value on the line
    y_on_line = m * XC + c;

    if YC > y_on_line
        position = 1;  % 'above line AB';
    elseif YC < y_on_line
        position = -1; % 'below line AB';
    else
        position = 0;  % 'on line AB';
    end
end

% % Example usage
% XA = 0; YA = 0;
% XB = 4; YB = 4;
% XC = 2; YC = 0;
% 
% [XD, YD, position] = find_perpendicular_point(XA, YA, XB, YB, XC, YC);
% fprintf('Point D is at (%.2f, %.2f)\n', XD, YD);
% fprintf('Point C is %s.\n', position);
