__precompile__()

module Master

using Types
export masternodes, gaussnodes, mkshape, mkmaster;

include("masternodes.jl");
include("gaussnodes.jl");
include("legendrepolynomial.jl");
include("simplexmonomial.jl");
include("tensorproduct.jl");
include("mkshape.jl");

struct MASTERStruct
    nd::IntP; # physical dimension
    porder::IntP; # polynomial degree
    pgauss::IntP; # order of Gauss quadrature
    elemtype::IntP; # type of elements
    nodetype::IntP; # type of nodes
    npe::IntP; # number of polynomials per element
    npf::IntP; # number of polynomials per face
    nge::IntP; # number of Gauss points per element
    ngf::IntP; # number of Gauss points per face
    np1d::IntP; # number of polynomials per 1D element
    ng1d::IntP; # number of Gauss points per 1D element
    perm::Array{IntP,2}; # to get face nodes from element nodes
    shapegt::Array{FloatP,3}; # element shape functions at Gauss points (transpose)
    shapegw::Array{FloatP,3}; # element shape functions at Gauss points multiplied by Gauss weights
    shapfgt::Array{FloatP,3}; # face shape functions at Gauss points (transpose)
    shapfgw::Array{FloatP,3}; # face shape functions at Gauss points multiplied by Gauss weights
    shapent::Array{FloatP,3}; # element shape functions at nodes (transpose)
    shapen::Array{FloatP,3};  # element shape functions at nodes
    shapfnt::Array{FloatP,3}; # face shape functions at nodes (transpose)
    shapfn::Array{FloatP,3};  # face shape functions at nodes
    xpe::Array{FloatP,2}; # nodal points on master element
    gpe::Array{FloatP,2}; # gauss points on master element
    gwe::Array{FloatP,2}; # gauss weighs on master element
    xpf::Array{FloatP,2}; # nodal points on master face
    gpf::Array{FloatP,2}; # gauss points on master face
    gwf::Array{FloatP,2}; # gauss weighs on master face
    shap1dgt::Array{FloatP,3}; # 1D shape functions at Gauss points (transpose)
    shap1dgw::Array{FloatP,3}; # 1D shape functions at Gauss points multiplied by Gauss weights
    shap1dnt::Array{FloatP,3}; # 1D shape functions at nodes (transpose)
    shap1dn::Array{FloatP,3}; # 1D shape functions at nodes
    xp1d::Array{FloatP,2}; # node points on 1D element
    gp1d::Array{FloatP,2}; # gauss points on 1D element
    gw1d::Array{FloatP,2}; # gauss weights on 1D element
    telem::Array{IntP,2}; # element
    tface::Array{IntP,2}; # face
end

function mkmaster(dim::IntP,porder::IntP,pgauss::IntP,elemtype::IntP,nodetype::IntP)

# node positions on the master element and master face
xpe,telem,xpf,tface,perm = masternodes(porder,dim,elemtype);

# Gauss points and weights on the master volume element
gpe, gwe = gaussnodes(pgauss,dim,elemtype);

# shape functions and derivatives on the master volume element
shapeg = mkshape(porder,xpe,gpe,elemtype);
shapen = mkshape(porder,xpe,xpe,elemtype);

npe,nge,nd1 = size(shapeg);
shapegt = zeros(nge,npe,nd1);
shapent = zeros(npe,npe,nd1);
shapegw = zeros(npe,nge,nd1);
Mgwe = zeros(nge,nge);
for i = 1:nge
    Mgwe[i,i] = gwe[i];
end
for d=1:nd1
    shapegt[:,:,d] = shapeg[:,:,d]';
    shapent[:,:,d] = shapen[:,:,d]';
    shapegw[:,:,d] = shapeg[:,:,d]*Mgwe;
end

if dim>1
    # Gauss points and weights on the master face element
    gpf, gwf = gaussnodes(pgauss,dim-1,elemtype);
    # shape functions and derivatives on the master face element
    shapfg = mkshape(porder,xpf,gpf,elemtype);
    shapfn = mkshape(porder,xpf,xpf,elemtype);
else
    gpf   = [0.0];
    gwf   = [1.0];
    shapfg = reshape([1.0],1,1,1);
    shapfn = reshape([1.0],1,1,1);
end

npf = size(shapfg,1);
ngf = size(shapfg,2);
shapfgt = zeros(ngf,npf,dim);
shapfnt = zeros(npf,npf,dim);
shapfgw = zeros(npf,ngf,dim);
Mgwf = zeros(ngf,ngf);
for i = 1:ngf
    Mgwf[i,i] = gwf[i];
end
for d=1:dim
    shapfgt[:,:,d] = shapfg[:,:,d]';
    shapfnt[:,:,d] = shapfn[:,:,d]';
    shapfgw[:,:,d] = shapfg[:,:,d]*Mgwf;
end

xp1d,~,~,~,~ = masternodes(porder,1,elemtype);
gp1d,gw1d = gaussnodes(pgauss,1,elemtype);
shap1dg = mkshape(porder,xp1d,gp1d,elemtype);
shap1dn = mkshape(porder,xp1d,xp1d,elemtype);

np1d = size(shap1dg,1);
ng1d = size(shap1dg,2);
shap1dgt = zeros(ng1d,np1d,2);
shap1dnt = zeros(np1d,np1d,2);
shap1dgw = zeros(np1d,ng1d,2);
Mgw1d = zeros(ng1d,ng1d);
for i = 1:ng1d
    Mgw1d[i,i] = gw1d[i];
end
for d=1:2
    shap1dgt[:,:,d] = shap1dg[:,:,d]';
    shap1dgw[:,:,d] = shap1dg[:,:,d]*Mgw1d;
    shap1dnt[:,:,d] = shap1dn[:,:,d]';
end

master = MASTERStruct(dim, porder, pgauss, elemtype, nodetype, npe, npf,
    nge, ngf, np1d, ng1d, perm, shapegt, shapegw, shapfgt, shapfgw, shapent, shapen,
    shapfnt, shapfn, xpe, gpe, gwe, xpf, gpf, gwf, shap1dgt, shap1dgw, shap1dnt,
    shap1dn, xp1d, gp1d, gw1d, telem, tface);

return master;

end

function mkmaster(app)
    master = mkmaster(app.nd,app.porder,app.pgauss,app.elemtype,app.nodetype);
    return master;
end

end
