function zb = localbathymetry(x, mu)
L = mu(2);
zbeach = mu(3);
zdeep = mu(4);

slope = (zdeep - zbeach)/L;
zb = zbeach + slope*x;
end
