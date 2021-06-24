
function [y,dydx] = meanLineNACAseries6(a,cl,x)

if nargin < 3
    x = linspace(0,1,100);
end

if a == 1; a = 0.9999; end

xEnd = x(end);
x = x/xEnd;

g = -(1/(1-a)) * (a^2*(0.5*log(a)-0.25) + 0.25);
h = (1/(1-a)) * (0.5*(1-a)^2*log(1-a) - 0.25*(1-a)^2) + g;

y = (cl/(2*pi*(a+1))) * ( (1/(1-a))* (0.5*(a-x).^2.*log(abs(a-x)) - 0.5*(1-x).^2.*log(1-x) + ...
        0.25*(1-x).^2 - 0.25*(a-x).^2) - x.*log(x) + g - h*x) * xEnd;

y(1) = 0;
y(end) = 0;

dydx = (cl/(2*pi*(a+1))) * ( (1/(1-a))* ((1-x).*log(1-x) - (a-x).*log(a-x) + 0.5*((a-x)-abs(a-x))) - 1 - log(x) - h);

% Correct last coordiate to avoid NaN if necessary:
if x(end) == 1
    eps = 1.e-3;
    dydx(end) = (cl/(2*pi*(a+1))) * ( (1/(1-a))* ((1-(x(end)-eps)).*log(1-(x(end)-eps)) - (a-(x(end)-eps)).*log(a-(x(end)-eps)) + 0.5*((a-(x(end)-eps))-abs(a-(x(end)-eps)))) - 1 - log((x(end)-eps)) - h);
end

end
