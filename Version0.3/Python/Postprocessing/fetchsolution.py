from numpy import *
from getsolution import getsolution

def fetchsolution(app,master,dmd,dirname="dataout"):

    nt = len(app['dt']);
    if (nt==1):
        Uout = getsolution(dirname + "/out",dmd,master['npe']);
    else:
        if len(app['soltime'])==0:
            app['soltime'] = arange(1,nt+1);
        fn = dirname + "/out_t" + str(app['soltime'][0]);
        tmp = getsolution(fn,dmd,master['npe']);
        if app['wave']==1:
            w = getsolution(dirname + "/out_wdg_t" + str(app['soltime'][0]),dmd,master['npe']);
            tmp = concatenate((tmp,w),axis=1);
        Uout = zeros((tmp.shape[0], tmp.shape[1], tmp.shape[2], len(app['soltime'])));
        Uout[:,:,:,0] = tmp;
        for i in range(1,len(app['soltime'])):
            fn = dirname + "/out_t" + str(app['soltime'][i]);
            tmp = getsolution(fn,dmd,master['npe']);
            if app['wave']==1:
                w = getsolution(dirname + "/out_wdg_t" + str(app['soltime'][i]),dmd,master['npe']);
                tmp = concatenate((tmp,w),axis=1);
            Uout[:,:,:,i] = tmp;
    return Uout
