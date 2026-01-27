function [fsw, dfsw] = switchfunctions(T, Tsw1, Tsw2, alpha)

sw1 = tanh(-alpha * (T - Tsw1) / pi ) * 0.5 + 0.5;
dsw1 = (alpha*(tanh((alpha*(T - Tsw1))/pi)^2 - 1))/(2*pi);

sw3 = tanh( alpha * (T - Tsw2) / pi ) * 0.5 + 0.5; 
dsw3 = -(alpha*(tanh((alpha*(T - Tsw2))/pi)^2 - 1))/(2*pi);

sw2 = 1.0 - sw1 - sw3;
dsw2 = -dsw1 - dsw3;

fsw = [sw1 sw2 sw3];
dfsw = [dsw1 dsw2 dsw3];
 
