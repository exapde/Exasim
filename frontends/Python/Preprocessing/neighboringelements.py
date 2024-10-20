from numpy import *

def neighboringelements(t2t, elem):

    t2te = t2t[:,elem];
    nbelem = sort(unique(t2te.flatten('F')));
    if nbelem[0] == -1:
        nbelem = nbelem[1:];

    return nbelem
