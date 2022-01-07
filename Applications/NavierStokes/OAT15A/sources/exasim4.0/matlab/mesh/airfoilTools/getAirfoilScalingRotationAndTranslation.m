
function [chord,theta,xTranslation,yTranslation] = getAirfoilScalingRotationAndTranslation()

fileName = 'scalingRotationAndTranslation_tmp.dat';

info = dlmread(fileName);

chord = info(1);
theta = info(2);
xTranslation = info(3);
yTranslation = info(4);

end
