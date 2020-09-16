from numpy import *

def getelemface(dim,elemtype):

    if dim==1:
        face = reshape([1,2],(2,1))
        nfe = 2
        nvf = 1
    elif dim==2:
        if elemtype==0:
            face = array([[2,3],[3,1],[1,2]])
            nfe = 3
            nvf = 2
        elif elemtype==1:
            face = array([[1,2],[2,3],[3,4],[4,1]])
            nfe = 4;
            nvf = 2;
    elif dim==3:
        if elemtype==0:
            face = array([[2,3,4],[1,4,3],[1,2,4],[1,3,2]])
            nfe = 4;
            nvf = 3;
        elif elemtype==1:
            face = array([[1,4,3,2],[5,6,7,8],[1,2,6,5],[3,4,8,7],[2,3,7,6],[4,1,5,8]])
            nfe = 6;
            nvf = 4;

    face = face.transpose() - 1;
    
    return face
