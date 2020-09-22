from numpy import *
from masternodes import *
from gaussnodes import *
# import importlib
# import mkshape #import the module here, so that it can be reloaded.
# importlib.reload(mkshape)
from mkshape import *

def mkmaster(dim,porder,pgauss,elemtype,nodetype):

    xpe,telem,xpf,tface,perm = masternodes(porder,dim,elemtype)

    gpe, gwe = gaussnodes(pgauss,dim,elemtype)
    gpe = array(gpe,float)
    gwe = array(gwe,float)
    gwe = reshape(gwe,(gwe.size,1))
    # shape functions and derivatives on the master volume element
    shapeg = mkshape(porder,xpe,gpe,elemtype)
    shapen = mkshape(porder,xpe,xpe,elemtype)

    npe     = shapeg.shape[0];
    nge     = shapeg.shape[1];
    nd1     = shapeg.shape[2];
    shapegt = zeros((nge,npe,nd1))
    shapent = zeros((npe,npe,nd1))
    shapegw = zeros((npe,nge,nd1))
    Mgwe = zeros((nge,nge))
    for i in range(0,nge):
        Mgwe[i,i] = gwe[i]

    for d in range(0,nd1):
        shapegt[:,:,d] = shapeg[:,:,d].transpose()
        shapent[:,:,d] = shapen[:,:,d].transpose()
        shapegw[:,:,d] = dot(shapeg[:,:,d],Mgwe)

    if dim>1:
        # Gauss points and weights on the master face element
        gpf, gwf = gaussnodes(pgauss,dim-1,elemtype)
        # shape functions and derivatives on the master face element
        shapfg = mkshape(porder,xpf,gpf,elemtype)
        shapfn = mkshape(porder,xpf,xpf,elemtype)
    else:
        gpf   = [0.0]
        gwf   = [1.0]
        shapfg =[1.0]
        shapfn =[1.0]

    npf = shapfg.shape[0];
    ngf = shapfg.shape[1];
    shapfgt = zeros((ngf,npf,dim));
    shapfnt = zeros((npf,npf,dim));
    shapfgw = zeros((npf,ngf,dim));
    Mgwf = zeros((ngf,ngf));
    for i in range(0,ngf):
        Mgwf[i,i] = gwf[i];
    for d in range(0,dim):
        shapfgt[:,:,d] = shapfg[:,:,d].transpose();
        shapfnt[:,:,d] = shapfn[:,:,d].transpose();
        shapfgw[:,:,d] = dot(shapfg[:,:,d],Mgwf);

    xp1d = masternodes(porder,1,elemtype)[0]
    gp1d,gw1d = gaussnodes(pgauss,1,elemtype)
    shap1dg = mkshape(porder,xp1d,gp1d,elemtype)
    shap1dn = mkshape(porder,xp1d,xp1d,elemtype)

    np1d = shap1dg.shape[0]
    ng1d = shap1dg.shape[1]
    shap1dgt = zeros((ng1d,np1d,2))
    shap1dnt = zeros((np1d,np1d,2))
    shap1dgw = zeros((np1d,ng1d,2))
    Mgw1d = zeros((ng1d,ng1d))
    for i in range(0,ng1d):
        Mgw1d[i,i] = gw1d[i]
    for d in range(0,2):
        shap1dgt[:,:,d] = shap1dg[:,:,d].transpose()
        shap1dgw[:,:,d] = dot(shap1dg[:,:,d],Mgw1d)
        shap1dnt[:,:,d] = shap1dn[:,:,d].transpose()

    master = {'nd' : dim, 'porder' : porder, 'pgauss' : pgauss, 'elemtype' : elemtype,
              'nodetype' : nodetype, 'npe' : npe, 'npf' : npf, 'nge' : nge, 'ngf' : ngf,
              'shapegt' : shapegt, 'shapegw' : shapegw, 'shapfgt' : shapfgt, 'shapfgw' : shapfgw,
              'shapent' : shapent, 'shapen' : shapen, 'shapfnt' : shapfnt, 'shapfn' : shapfn,
              'xpe' : xpe, 'gpe' : gpe, 'gwe' : gwe, 'xpf' : xpf, 'gpf' : gpf, 'gwf' : gwf,
              'shap1dgt' : shap1dgt, 'shap1dgw' : shap1dgw, 'shap1dnt' : shap1dnt,
              'shap1dn' : shap1dn, 'xp1d' : xp1d, 'gp1d' : gp1d, 'gw1d' : gw1d,
              'telem' : telem, 'tface' : tface, 'perm' : perm};

    return master

    # master = {'dim' : [], 'porder' : [], 'pgauss' : [], 'elemtype' : [], 'nodetype' : [],
    #           'npe' : [], 'npf' : [], 'nge' : [], 'ngf' : [], 'shapegt' : [], 'shapegw' : [],
    #           'shapfgt' : [], 'shapfgw' : [], 'shapent' : [], 'shapen' : [], 'shapfnt' : [],
    #           'shapfn' : [], 'xpe' : [], 'gpe' : [], 'gwe' : [], 'xpf' : [], 'gpf' : [],
    #           'gwf' : [], 'shap1dgt' : [], 'shap1dgw' : [], 'shap1dnt' : [], 'shap1dn' : [],
    #           'xp1d' : [], 'gp1d' : [], 'gw1d' : [], 'telem' : [], 'tface' : []};
    #
    # master['dim'] = dim
    # master['porder'] = porder
    # master['pgauss'] = pgauss
    # master['elemtype'] = elemtype
    # master['nodetype'] = nodetype
    # master['xpe'] = xpe
    # master['telem'] = telem
    # master['xpf'] = xpf
    # master['tface'] = tface
    # master['perm'] = perm
    # master['shapen'] = shapen;
    # master['xp1d'] = xp1d
    # shapvl = mkshape(porder,plocal,gpvl,elemtype)
    # master['shapvl'] = shapvl
    #
    # if ndim>1:
    #     gpfc,gwfc = gaussquad(pgauss,ndim-1,elemtype)
    #     master['gpfc'] = gpfc
    #     master['gwfc'] = gwfc
    #
    #     shapfc = mkshape(porder,plocfc,gpfc,elemtype)
    #     master['shapfc'] = shapfc
    #     gpfc = array(gpfc,float)
    #     gwfc = array(gwfc,float)
    #     master['gpfc'] = gpfc
    #     master['gwfc'] = reshape(gwfc,(gwfc.size,1))
    # else:
    #     master['gpfc'] = 0.5
    #     master['gwfc'] = 1.0
    #     master['shapfc'] = 1.0
    #
    # npv = master['shapvl'].shape[0]
    # ngv = master['shapvl'].shape[1]
    # master['shapvgdotshapvl']  = zeros((npv*npv,ngv,ndim+1))
    # master['shapvt'] = zeros((ngv,npv,ndim+1))
    # master['shapvg'] = zeros((npv,ngv,ndim+1))
    # for d in range(0,ndim+1):
    #     master['shapvt'][:,:,d] = master['shapvl'][:,:,d].transpose()
    #     master['shapvg'][:,:,d] = dot(master['shapvl'][:,:,d],diag(master['gwvl'][:,0]))
    #     for ii in range(0,npv):
    #         for jj in range(0,npv):
    #             master['shapvgdotshapvl'][ii*npv+jj,:,d] = master['shapvg'][jj,:,d]*master['shapvl'][ii,:,0]
    #
    #
    # npf = master['shapfc'].shape[0]
    # ngf = master['shapfc'].shape[1]
    # master['shapfgdotshapfc']  = zeros((npf*npf,ngf,ndim))
    # master['shapft'] = zeros((ngf,npf,ndim))
    # master['shapfg'] = zeros((npf,ngf,ndim))
    # for d in range(0,ndim):
    #     master['shapft'][:,:,d] = master['shapfc'][:,:,d].transpose()
    #     master['shapfg'][:,:,d] = dot(master['shapfc'][:,:,d],diag(master['gwfc'][:,0]))
    #     for ii in range(0,npf):
    #         for jj in range(0,npf):
    #             master['shapfgdotshapfc'][ii*npf+jj,:,d] = master['shapfg'][jj,:,d]*master['shapfc'][ii,:,0]
    #
    # master['shapmv'] = master['shapvt']
    # master['shapmf'] = master['shapft']
    # master['shapmh'] = master['shapmf']
    # master['permgeom'] = master['perm']


# function mkmaster(dim::IntP,porder::IntP,pgauss::IntP,elemtype::IntP,nodetype::IntP)
#
# # node positions on the master element and master face
# xpe,telem,xpf,tface,perm = masternodes(porder,dim,elemtype);
#
# # Gauss points and weights on the master volume element
# gpe, gwe = gaussnodes(pgauss,dim,elemtype);
#
# # shape functions and derivatives on the master volume element
# shapeg = mkshape(porder,xpe,gpe,elemtype);
# shapen = mkshape(porder,xpe,xpe,elemtype);
#
# npe,nge,nd1 = size(shapeg);
# shapegt = zeros(nge,npe,nd1);
# shapent = zeros(npe,npe,nd1);
# shapegw = zeros(npe,nge,nd1);
# Mgwe = zeros(nge,nge);
# for i = 1:nge
#     Mgwe[i,i] = gwe[i];
# end
# for d=1:nd1
#     shapegt[:,:,d] = shapeg[:,:,d]';
#     shapent[:,:,d] = shapen[:,:,d]';
#     shapegw[:,:,d] = shapeg[:,:,d]*Mgwe;
# end
#
# if dim>1
#     # Gauss points and weights on the master face element
#     gpf, gwf = gaussnodes(pgauss,dim-1,elemtype);
#     # shape functions and derivatives on the master face element
#     shapfg = mkshape(porder,xpf,gpf,elemtype);
#     shapfn = mkshape(porder,xpf,xpf,elemtype);
# else
#     gpf   = [0.0];
#     gwf   = [1.0];
#     shapfg = reshape([1.0],1,1,1);
#     shapfn = reshape([1.0],1,1,1);
# end
#
# npf = size(shapfg,1);
# ngf = size(shapfg,2);
# shapfgt = zeros(ngf,npf,dim);
# shapfnt = zeros(npf,npf,dim);
# shapfgw = zeros(npf,ngf,dim);
# Mgwf = zeros(ngf,ngf);
# for i = 1:ngf
#     Mgwf[i,i] = gwf[i];
# end
# for d=1:dim
#     shapfgt[:,:,d] = shapfg[:,:,d]';
#     shapfnt[:,:,d] = shapfn[:,:,d]';
#     shapfgw[:,:,d] = shapfg[:,:,d]*Mgwf;
# end
#
# xp1d,~,~,~,~ = masternodes(porder,1,elemtype);
# gp1d,gw1d = gaussnodes(pgauss,1,elemtype);
# shap1dg = mkshape(porder,xp1d,gp1d,elemtype);
# shap1dn = mkshape(porder,xp1d,xp1d,elemtype);
#
# np1d = size(shap1dg,1);
# ng1d = size(shap1dg,2);
# shap1dgt = zeros(ng1d,np1d,2);
# shap1dnt = zeros(np1d,np1d,2);
# shap1dgw = zeros(np1d,ng1d,2);
# Mgw1d = zeros(ng1d,ng1d);
# for i = 1:ng1d
#     Mgw1d[i,i] = gw1d[i];
# end
# for d=1:2
#     shap1dgt[:,:,d] = shap1dg[:,:,d]';
#     shap1dgw[:,:,d] = shap1dg[:,:,d]*Mgw1d;
#     shap1dnt[:,:,d] = shap1dn[:,:,d]';
# end
#
# master = Structs.MasterStruct(dim, porder, pgauss, elemtype, nodetype, npe, npf,
#     nge, ngf, np1d, ng1d, perm, shapegt, shapegw, shapfgt, shapfgw, shapent, shapen,
#     shapfnt, shapfn, xpe, gpe, gwe, xpf, gpf, gwf, shap1dgt, shap1dgw, shap1dnt,
#     shap1dn, xp1d, gp1d, gw1d);
#
# return master;
#
# end
#
# function mkmaster(app)
#     master = mkmaster(app.nd,app.porder,app.pgauss,app.elemtype,app.nodetype);
#     return master;
# end
#
# end
