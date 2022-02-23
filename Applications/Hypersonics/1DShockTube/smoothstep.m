function [out] = smoothstep(x,y1,y2,b)
if y1 < y2
    flip = 1;
    ymin = y1;
    ymax = y2;
else
    flip = -1;
    ymin = y2;
    ymax = y1;
end
height = ymax - ymin;
scale = (height/2);
out = flip*scale*tanh(x/b);
disp(scale + min(out(:)));
% out = out + (ymin - min(out(:)));
out = out + (ymin + scale);
end
