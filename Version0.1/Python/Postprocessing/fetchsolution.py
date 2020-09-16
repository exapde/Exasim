from numpy import *
from getsolution import getsolution

def fetchsolution(app,master,dmd):

    nt = len(app['dt']);
    if (nt==1):
        Uout = getsolution("dataout/out",dmd,master['npe']);
    else:
        if len(app['vistime'])==0:
            app['vistime'] = arange(1,nt);
        fn = "dataout/out_t" + str(app['vistime'][0]);
        tmp = getsolution(fn,dmd,master['npe']);
        Uout = zeros((tmp.shape[0], tmp.shape[1], tmp.shape[2], len(app['vistime'])));
        Uout[:,:,:,0] = tmp;
        tmp = [];
        for i in range(1,len(app['vistime'])):
            fn = "dataout/out_t" + str(app['vistime'][i]);
            Uout[:,:,:,i] = getsolution(fn,dmd,master['npe']);
    return Uout
