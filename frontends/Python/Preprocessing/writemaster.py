from numpy import *

def writemaster(master,filename):

    ndims = zeros((20,1));
    ndims[1-1] = master['nd'];
    ndims[2-1] = master['elemtype'];
    ndims[3-1] = master['nodetype'];
    ndims[4-1] = master['porder'];
    ndims[5-1] = master['pgauss'];
    ndims[6-1] = master['npe'];
    ndims[7-1] = master['npf'];
    ndims[8-1] = master['nge'];
    ndims[9-1] = master['ngf'];
    ndims[10-1] = size(master['xp1d']);
    ndims[11-1] = size(master['gp1d']);

    nsize = zeros((22,1));
    nsize[1-1] = size(ndims);
    nsize[2-1] = size(master['shapegt']);
    nsize[3-1] = size(master['shapegw']);
    nsize[4-1] = size(master['shapfgt']);
    nsize[5-1] = size(master['shapfgw']);
    nsize[6-1] = size(master['shapent']);
    nsize[7-1] = size(master['shapen']);
    nsize[8-1] = size(master['shapfnt']);
    nsize[9-1] = size(master['shapfn']);
    nsize[10-1] = size(master['xpe']);
    nsize[11-1] = size(master['gpe']);
    nsize[12-1] = size(master['gwe']);
    nsize[13-1] = size(master['xpf']);
    nsize[14-1] = size(master['gpf']);
    nsize[15-1] = size(master['gwf']);
    nsize[16-1] = size(master['shap1dgt']);
    nsize[17-1] = size(master['shap1dgw']);
    nsize[18-1] = size(master['shap1dnt']);
    nsize[19-1] = size(master['shap1dn']);
    nsize[20-1] = size(master['xp1d']);
    nsize[21-1] = size(master['gp1d']);
    nsize[22-1] = size(master['gw1d']);

    print("Writing master into file...");
    fileID = open(filename, 'wb');

    # write master structure to files
    array(size(nsize), dtype=float64).tofile(fileID)
    nsize.astype('float64').tofile(fileID)
    ndims.astype('float64').tofile(fileID);
    master['shapegt'] = array(master['shapegt']).flatten(order = 'F');
    master['shapegt'].astype('float64').tofile(fileID);
    master['shapegw'] = array(master['shapegw']).flatten(order = 'F');
    master['shapegw'].astype('float64').tofile(fileID);
    master['shapfgt'] = array(master['shapfgt']).flatten(order = 'F');
    master['shapfgt'].astype('float64').tofile(fileID);
    master['shapfgw'] = array(master['shapfgw']).flatten(order = 'F');
    master['shapfgw'].astype('float64').tofile(fileID);
    master['shapent'] = array(master['shapent']).flatten(order = 'F');
    master['shapent'].astype('float64').tofile(fileID);
    master['shapen'] = array(master['shapen']).flatten(order = 'F');
    master['shapen'].astype('float64').tofile(fileID);
    master['shapfnt'] = array(master['shapfnt']).flatten(order = 'F');
    master['shapfnt'].astype('float64').tofile(fileID);
    master['shapfn'] = array(master['shapfn']).flatten(order = 'F');
    master['shapfn'].astype('float64').tofile(fileID);
    master['xpe'] = array(master['xpe']).flatten(order = 'F');
    master['xpe'].astype('float64').tofile(fileID);
    master['gpe'] = array(master['gpe']).flatten(order = 'F');
    master['gpe'].astype('float64').tofile(fileID);
    master['gwe'] = array(master['gwe']).flatten(order = 'F');
    master['gwe'].astype('float64').tofile(fileID);
    master['xpf'] = array(master['xpf']).flatten(order = 'F');
    master['xpf'].astype('float64').tofile(fileID);
    master['gpf'] = array(master['gpf']).flatten(order = 'F');
    master['gpf'].astype('float64').tofile(fileID);
    master['gwf'] = array(master['gwf']).flatten(order = 'F');
    master['gwf'].astype('float64').tofile(fileID);
    master['shap1dgt'] = array(master['shap1dgt']).flatten(order = 'F');
    master['shap1dgt'].astype('float64').tofile(fileID);
    master['shap1dgw'] = array(master['shap1dgw']).flatten(order = 'F');
    master['shap1dgw'].astype('float64').tofile(fileID);
    master['shap1dnt'] = array(master['shap1dnt']).flatten(order = 'F');
    master['shap1dnt'].astype('float64').tofile(fileID);
    master['shap1dn'] = array(master['shap1dn']).flatten(order = 'F');
    master['shap1dn'].astype('float64').tofile(fileID);
    master['xp1d'] = array(master['xp1d']).flatten(order = 'F');
    master['xp1d'].astype('float64').tofile(fileID);
    master['gp1d'] = array(master['gp1d']).flatten(order = 'F');
    master['gp1d'].astype('float64').tofile(fileID);
    master['gw1d'] = array(master['gw1d']).flatten(order = 'F');
    master['gw1d'].astype('float64').tofile(fileID);

    fileID.close();
