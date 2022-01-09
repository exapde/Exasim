from numpy import *
from localbasis import localbasis
from getelemface import getelemface
from permindex import permindex
import sys

def faceconnectivity2(t,f2t,dim,elemtype,porder):

    nf = f2t.shape[1];
    philocvl,philocfc,tm1,plocfc,perm = localbasis(porder,dim,elemtype)[0:5]
    perm = perm-1;
    permind = permindex(plocfc,dim,elemtype);
    face = getelemface(dim,elemtype);

    npf = perm.shape[0];
    npe = philocvl.shape[0];

    facecon = zeros((npf,2,nf)).astype(int);
    if dim<=2:
        for i in range(0,nf):
            e1 = f2t[0,i]-1;
            l1 = f2t[1,i]-1;
            e2 = f2t[2,i]-1;
            l2 = f2t[3,i]-1;

            facecon[:,0,i] = e1*npe + perm[:,l1];

            if e2>=0: # face i is an interior face
                facecon[:,1,i] = e2*npe + perm[permind,l2];
            else:
                facecon[:,1,i] = facecon[:,0,i];
    else:
        for i in range(0,nf):
            e1 = f2t[0,i]-1;
            l1 = f2t[1,i]-1;
            e2 = f2t[2,i]-1;
            l2 = f2t[3,i]-1;

            facecon[:,0,i] = e1*npe + perm[:,l1];

            if e2>=0: # face i is an interior face
                f1 = t[face[:,l1],e1];
                f2 = t[face[:,l2],e2];
                if elemtype==0:
                    if (f1[0]==f2[0]) and (f1[1]==f2[2]) and (f1[2]==f2[1]):
                        k = 0;
                    elif (f1[0]==f2[1]) and (f1[1]==f2[0]) and (f1[2]==f2[2]):
                        k = 1;
                    elif (f1[0]==f2[2]) and (f1[1]==f2[1]) and (f1[2]==f2[0]):
                        k = 2;
                    else:
                        error("Mesh connectivity is wrong");
                else:
                    if (f1[0]==f2[0]) and (f1[1]==f2[3]) and (f1[2]==f2[2]) and (f1[3]==f2[1]):
                        k = 0;
                    elif (f1[0]==f2[3]) and (f1[1]==f2[2]) and (f1[2]==f2[1]) and (f1[3]==f2[0]):
                        k = 1;
                    elif (f1[0]==f2[2]) and (f1[1]==f2[1]) and (f1[2]==f2[0]) and (f1[3]==f2[3]):
                        k = 2;
                    elif (f1[0]==f2[1]) and (f1[1]==f2[0]) and (f1[2]==f2[3]) and (f1[3]==f2[2]):
                        k = 3;
                    else:
                        error("Mesh connectivity is wrong");

                facecon[:,1,i] = e2*npe + perm[permind[:,k],l2];
            else:
                facecon[:,1,i] = facecon[:,0,i];

    return facecon
