function xshift = shifted_curve_normal_direction(x, d)

n = size(x,1);
xshift = 0*x;
for i = 1:n-2
  XA = x(i+1,1);
  YA = x(i+1,2);
  s1 = (YA - x(i,2)) / (XA - x(i,1));
  s2 = (x(i+2,2) - YA) / (x(i+2,1) - XA);
  s = 0.5*(s1+s2);
  p = -1 / s;  
  
 % Compute the offsets
  dx = sign(s) * d / sqrt(1 + p^2);
  dy = p * dx;  
  
  xshift(i+1,1) = XA + dx;
  xshift(i+1,2) = YA + dy;
end

i = 1;
XA = x(i+1,1);
YA = x(i+1,2);
s = (YA - x(i,2)) / (XA - x(i,1));
p = -1 / s;
dx = d / sqrt(1 + p^2);
dy = p * dx;  
xshift(i,1) = XA + dx;
xshift(i,2) = YA + dy;

i = n;
XA = x(i,1);
YA = x(i,2);
s = (YA - x(i-1,2)) / (XA - x(i-1,1));
p = -1 / s;
dx = d / sqrt(1 + p^2);
dy = p * dx;  
xshift(i,1) = XA + dx;
xshift(i,2) = YA + dy;

