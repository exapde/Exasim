function y=loginc(x,alpha)
a=min(x(:));b=max(x(:));
y = a + (b-a)*(exp(alpha*(x-a)/(b-a))-1)/(exp(alpha)-1);
