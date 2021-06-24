
function dist = dist_NACA65_X10_afterScalingAndRotation(p,X,a)
% Distance function for NACA 65-X10 airfoil, with scaling, rotation
% and translation applied to form the rotor blade for the PW problem.
%

% chord = 1;  % chord length scaling
% theta = 60.8*pi/180; % anticlockwise rotation angle applied to airfoil
% TRANS = 0;  % x-translation applied to scaled & rotated sock mesh

if nargin<3; a = 1.0; end

[chord,theta,xTranslation,yTranslation] = getAirfoilScalingRotationAndTranslation();

[xf,yf] = NACA65_X10(X,a);

% Scale
pf = chord*[xf yf];

% Rotate (counterclockwise)
A_ROT = [cos(theta) -sin(theta); sin(theta) cos(theta)];
pf = (A_ROT*(pf'))';

% Translate
pf(:,1) = pf(:,1) + xTranslation;
pf(:,2) = pf(:,2) + yTranslation;

% Distance
dist = dpoly(p,pf);

end
