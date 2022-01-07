
function writeAirfoilScalingRotationAndTranslation(chord,theta,xTranslation,yTranslation)

fileName = 'scalingRotationAndTranslation_tmp.dat';

info = [chord,theta,xTranslation,yTranslation];

dlmwrite(fileName,info);

end
