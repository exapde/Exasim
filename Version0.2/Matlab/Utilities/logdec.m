function y=logdec(x,alpha)
a=min(x(:));b=max(x(:));
y = a + (b-a)*(1-exp(-alpha*(x-a)/(b-a)))/(1-exp(-alpha));
