
function [xf,yf] = NACA65_X10(X,a,minThickness)

if nargin < 3; minThickness = 0; end

if X~=18
    fileName = 'n65-010.dat';

    cl = X/10;      % Lift coefficient.

    coordinatesSymmetric = dlmread(fileName);

    numPoints = size(coordinatesSymmetric,1);

    xSymmetric = coordinatesSymmetric(:,1);
    xAlongChord = xSymmetric(ceil(numPoints/2)+1:end);

    thickness_tmp = coordinatesSymmetric(:,2);
    thickness = abs(thickness_tmp(ceil(numPoints/2)+1:end));

    coordToChange = find(thickness < minThickness);
    thickness(coordToChange) = minThickness;

    [yMean,dydx] = meanLineNACAseries6(a,cl,xAlongChord);

    % slope = [(yMean(2) - yMean(1))/(xAlongChord(2) - xAlongChord(1)); (yMean(3:end) - yMean(1:end-2))./(xAlongChord(3:end) - xAlongChord(1:end-2)); (yMean(end) - yMean(end-1))/(xAlongChord(end) - xAlongChord(end-1))];

    theta = atan(dydx);

    yBottom = yMean - thickness.*cos(theta);
    xBottom = xAlongChord + thickness.*sin(theta);

    yUpper = yMean + thickness.*cos(theta);
    xUpper = xAlongChord - thickness.*sin(theta);

    % yBottom = yBottom(2:end);
    % xBottom = xBottom(2:end);
    % 
    % yUpper = yUpper(2:end);
    % xUpper = xUpper(2:end);

    xf = [flip(xBottom); xUpper];
    yf = [flip(yBottom); yUpper];

    % If minThickness > 0, close airfoil surface.
    if minThickness > 0
        xTE = 0.5*(xf(1)+xf(end));
        yTE = 0.5*(yf(1)+yf(end));

        xf(1) = xTE;
        yf(1) = yTE;

        xf(end) = xTE;
        yf(end) = yTE;
    end

    % xf = [flip(xBottom); 0; xUpper];
    % yf = [flip(yBottom); 0; yUpper];
    
else
    [xf,yf] = read_foil('n651810');
end
