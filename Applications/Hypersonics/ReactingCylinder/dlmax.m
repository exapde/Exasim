function dout = dlmax(x, alpha)
    dout = atan(alpha*(x))/pi + (alpha*(x))./(pi*(alpha^2*(x).^2 + 1)) + 1/2;
end