from numpy import *
from xiny import xiny

def permindex(plocfc,dim,elemtype):

    npf = plocfc.shape[0];

    if dim==1:
        ind = 0;
    elif dim==2:
        ind = arange(npf-1,-1,-1);
    elif dim==3:
        if elemtype==0:
            ind = zeros((npf,3)).astype(int);

            # [1 3 2]
            plocfc2 = 0*plocfc;
            plocfc2[:,0] = plocfc[:,1];
            plocfc2[:,1] = plocfc[:,0];
            ind[:,0] = xiny(plocfc.round(8),plocfc2.round(8),1);

            # [2 1 3]
            plocfc2[:,0] = 1.0 - plocfc[:,0] - plocfc[:,1];
            plocfc2[:,1] = plocfc[:,1];
            ind[:,1] = xiny(plocfc.round(8),plocfc2.round(8),1);

            # [3 2 1]
            plocfc2[:,0] = plocfc[:,0];
            plocfc2[:,1] = 1.0 - plocfc[:,0] - plocfc[:,1];
            ind[:,2] = xiny(plocfc.round(8),plocfc2.round(8),1);

        # plocfc2 = plocfc;
        # plocfc2(:,1) = plocfc(:,2);
        # plocfc2(:,2) = plocfc(:,1);
        # ind(:,1) = xiny(round(plocfc,8),round(plocfc2,8));
        #
        # % [2 1 3]
        # plocfc2 = plocfc;
        # plocfc2(:,1) = 1-plocfc(:,1)-plocfc(:,2);
        # ind(:,2) = xiny(round(plocfc,8),round(plocfc2,8));
        #
        # % [3 2 1]
        # plocfc2 = plocfc;
        # plocfc2(:,2) = 1-plocfc(:,1)-plocfc(:,2);
        # ind(:,3) = xiny(round(plocfc,8),round(plocfc2,8));                                        
        else:
            ind = zeros((npf,4)).astype(int);

            # [1 4 3 2]
            plocfc2 = 0*plocfc;
            plocfc2[:,0] = plocfc[:,1];
            plocfc2[:,1] = plocfc[:,0];
            ind[:,0] = xiny(plocfc.round(8),plocfc2.round(8),1);

            # [2 1 4 3]
            plocfc2[:,0] = plocfc[:,0];
            plocfc2[:,1] = 1.0 - plocfc[:,1];
            ind[:,1] = xiny(plocfc.round(8),plocfc2.round(8),1);

            # [3 2 1 4]
            plocfc2[:,0] = 1.0 - plocfc[:,1];
            plocfc2[:,1] = 1.0 - plocfc[:,0];
            ind[:,2] = xiny(plocfc.round(8),plocfc2.round(8),1);

            # [4 3 2 1]
            plocfc2[:,0] = 1.0 - plocfc[:,0];
            plocfc2[:,1] = plocfc[:,1];
            ind[:,3] = xiny(plocfc.round(8),plocfc2.round(8),1);

    return ind
