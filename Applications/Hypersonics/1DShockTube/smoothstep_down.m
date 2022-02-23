function [out] = smoothstep_down(x,yminus,yplus,b)
% if yminus < yplus
%     flip = 1;
%     ymin = yminus;
%     ymax = yplus;
% else
    flip = -1;
    ymin = yplus;
    ymax = yminus;
% end
height = ymax - ymin;
scale = (height/2);
out = flip*scale*tanh(x/b);
% out = out + (ymin - min(out(:)));
out = out + (ymin + scale);
end
