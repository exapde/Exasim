from numpy import *

def writebin(filename,a):

    fileID = open(filename, 'wb');
    a.flatten(order = 'F').astype('float64').tofile(fileID);
    fileID.close();
    return 0;
