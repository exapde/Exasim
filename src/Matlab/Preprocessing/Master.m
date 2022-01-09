function master = Master(app)

dim = app.nd;
porder = app.porder;
pgauss = app.pgauss;
elemtype = app.elemtype;
nodetype = app.nodetype;

% node positions on the master element and master face
[xpe,telem,xpf,tface,perm] = masternodes(porder,dim,elemtype);
if (dim==1)
    xpf = 0;
    tface = 1;
end

% Gauss points and weights on the master volume element
[gpe, gwe] = gaussnodes(pgauss,dim,elemtype);

% shape functions and derivatives on the master volume element
shapeg = mkshape(porder,xpe,gpe,elemtype);
shapen = mkshape(porder,xpe,xpe,elemtype);

[npe,nge,nd1] = size(shapeg);
shapegt = zeros(nge,npe,nd1);
shapent = zeros(npe,npe,nd1);
shapegw = zeros(npe,nge,nd1);
Mgwe = zeros(nge,nge);
for i = 1:nge
    Mgwe(i,i) = gwe(i);
end
for d=1:nd1
    shapegt(:,:,d) = shapeg(:,:,d)';
    shapent(:,:,d) = shapen(:,:,d)';
    shapegw(:,:,d) = shapeg(:,:,d)*Mgwe;
end

if dim>1
    % Gauss points and weights on the master face element
    [gpf, gwf] = gaussnodes(pgauss,dim-1,elemtype);
    % shape functions and derivatives on the master face element
    shapfg = mkshape(porder,xpf,gpf,elemtype);
    shapfn = mkshape(porder,xpf,xpf,elemtype);
else
    gpf   = 0.0;
    gwf   = 1.0;
    shapfg = reshape(1.0,1,1,1);
    shapfn = reshape(1.0,1,1,1);
end

npf = size(shapfg,1);
ngf = size(shapfg,2);
shapfgt = zeros(ngf,npf,dim);
shapfnt = zeros(npf,npf,dim);
shapfgw = zeros(npf,ngf,dim);
Mgwf = zeros(ngf,ngf);
for i = 1:ngf
    Mgwf(i,i) = gwf(i);
end
for d=1:dim
    shapfgt(:,:,d) = shapfg(:,:,d)';
    shapfnt(:,:,d) = shapfn(:,:,d)';
    shapfgw(:,:,d) = shapfg(:,:,d)*Mgwf;
end

xp1d = masternodes(porder,1,elemtype);
[gp1d,gw1d] = gaussnodes(pgauss,1,elemtype);
shap1dg = mkshape(porder,xp1d,gp1d,elemtype);
shap1dn = mkshape(porder,xp1d,xp1d,elemtype);

np1d = size(shap1dg,1);
ng1d = size(shap1dg,2);
shap1dgt = zeros(ng1d,np1d,2);
shap1dnt = zeros(np1d,np1d,2);
shap1dgw = zeros(np1d,ng1d,2);
Mgw1d = zeros(ng1d,ng1d);
for i = 1:ng1d
    Mgw1d(i,i) = gw1d(i);
end
for d=1:2
    shap1dgt(:,:,d) = shap1dg(:,:,d)';
    shap1dgw(:,:,d) = shap1dg(:,:,d)*Mgw1d;
    shap1dnt(:,:,d) = shap1dn(:,:,d)';
end

master.dim = dim;
master.porder = porder;
master.pgauss = pgauss;
master.elemtype = elemtype;
master.nodetype = nodetype;
master.npe = npe;
master.npf = npf;
master.nge = nge;
master.ngf = ngf;
master.np1d = np1d;
master.ng1d = ng1d;
master.perm = perm;
master.shapegt = shapegt;
master.shapegw = shapegw;
master.shapfgt = shapfgt;
master.shapfgw = shapfgw;
master.shapent = shapent;
master.shapen = shapen;
master.shapfnt = shapfnt;
master.shapfn = shapfn;
master.xpe = xpe;
master.gpe = gpe;
master.gwe = gwe;
master.xpf = xpf;
master.gpf = gpf;
master.gwf = gwf;
master.shap1dgt = shap1dgt;
master.shap1dgw = shap1dgw;
master.shap1dnt = shap1dnt;
master.shap1dn = shap1dn;
master.xp1d = xp1d;
master.gp1d = gp1d;
master.gw1d = gw1d;
master.telem = telem;
master.tface = tface;


