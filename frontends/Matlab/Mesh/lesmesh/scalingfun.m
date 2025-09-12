function F = scalingfun(x,n,c)
     
F = c;
for i = 1:1:(n-1)
    F = F + x.^i;
end
