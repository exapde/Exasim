function localbasis(porder,dim,elemtype)

pelem,telem,pface,tface,perm = masternodes(porder,dim,elemtype);

npe = size(pelem,1);
npf = size(pface,1);

if dim==1
    nve = 2;
    nvf = 1;
elseif dim==2 && elemtype==0 # tri
    nve = 3;
    nvf = 2;
elseif dim==2 && elemtype==1 # quad
    nve = 4;
    nvf = 2;
elseif dim==3 && elemtype==0 # tet
    nve = 4;
    nvf = 3;
elseif dim==3 && elemtype==1 # hex
    nve = 8;
    nvf = 4;
end

phielem = zeros(Float64, npe, nve);
phiface = zeros(Float64, npf, nvf);

if dim==1
    xi  = pelem[:,1];
    phielem[:,1] = 1.0 .- xi;
    phielem[:,2] = xi;
    phiface = [1.0];
elseif dim==2 && elemtype==0 # tri
    xi  = pelem[:,1];
    eta = pelem[:,2];
    phielem[:,1] = 1.0 .- xi .- eta;
    phielem[:,2] = xi;
    phielem[:,3] = eta;

    xi = pface[:,1];
    phiface[:,1] = 1.0 .- xi;
    phiface[:,2] = xi;
elseif dim==2 && elemtype==1 # quad
    xi  = pelem[:,1];
    eta = pelem[:,2];
    phielem[:,1] = (1.0 .- xi).*(1.0 .- eta);
    phielem[:,2] = xi.*(1.0 .- eta);
    phielem[:,3] = xi.*eta;
    phielem[:,4] = (1.0 .- xi).*eta;

    xi = pface[:,1];
    phiface[:,1] = 1.0 .- xi;
    phiface[:,2] = xi;
elseif dim==3 && elemtype==0 # tet
    xi   = pelem[:,1];
    eta  = pelem[:,2];
    zeta = pelem[:,3];
    phielem[:,1] = 1.0 .- xi .- eta .- zeta;
    phielem[:,2] = xi;
    phielem[:,3] = eta;
    phielem[:,4] = zeta;

    xi = pface[:,1];
    eta = pface[:,2];
    phiface[:,1] = 1.0 .- xi .- eta;
    phiface[:,2] = xi;
    phiface[:,3] = eta;
elseif dim==3 && elemtype==1 # hex
    xi   = pelem[:,1];
    eta  = pelem[:,2];
    zeta = pelem[:,3];
    phielem[:,1] = (1.0 .-xi).*(1.0 .-eta).*(1.0 .-zeta);
    phielem[:,2] = xi.*(1.0 .-eta).*(1.0 .-zeta);
    phielem[:,3] = xi.*eta.*(1.0 .-zeta);
    phielem[:,4] = (1.0 .-xi).*eta.*(1.0 .-zeta);
    phielem[:,5] = (1.0 .-xi).*(1.0 .-eta).*(zeta);
    phielem[:,6] = xi.*(1.0 .-eta).*(zeta);
    phielem[:,7] = xi.*eta.*(zeta);
    phielem[:,8] = (1.0 .-xi).*eta.*(zeta);

    xi = pface[:,1];
    eta = pface[:,2];
    phiface[:,1] = (1.0 .- xi).*(1.0 .- eta);
    phiface[:,2] = xi.*(1.0 .- eta);
    phiface[:,3] = xi.*eta;
    phiface[:,4] = (1.0 .- xi).*eta;
end

return phielem,phiface,pelem,pface,perm

end
