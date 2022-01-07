from numpy import *

def mymax(a):
    if a.size==0:
        maxa=0;
    else:
        maxa = max(a);
    return maxa;


def checkmesh(mesh1,mesh2):

    diff = ones(22);
    diff[0] = mymax(abs(mesh1['nsize'].flatten('F')-mesh2['nsize'].flatten('F')));
    diff[1] = mymax(abs(mesh1['ndims'].flatten('F')-mesh2['ndims'].flatten('F')));
    diff[2] = mymax(abs(mesh1['facecon'].flatten('F')-mesh2['facecon'].flatten('F')));
    diff[3] = mymax(abs(mesh1['eblks'].flatten('F')-mesh2['eblks'].flatten('F')));
    diff[4] = mymax(abs(mesh1['fblks'].flatten('F')-mesh2['fblks'].flatten('F')));
    diff[5] = mymax(abs(mesh1['nbsd'].flatten('F')-mesh2['nbsd'].flatten('F')));
    diff[6] = mymax(abs(mesh1['elemsend'].flatten('F')-mesh2['elemsend'].flatten('F')));
    diff[7] = mymax(abs(mesh1['elemrecv'].flatten('F')-mesh2['elemrecv'].flatten('F')));
    diff[8] = mymax(abs(mesh1['elemsendpts'].flatten('F')-mesh2['elemsendpts'].flatten('F')));
    diff[9] = mymax(abs(mesh1['elemrecvpts'].flatten('F')-mesh2['elemrecvpts'].flatten('F')));
    diff[10] = mymax(abs(mesh1['elempart'].flatten('F')-mesh2['elempart'].flatten('F')));
    diff[11] = mymax(abs(mesh1['elempartpts'].flatten('F')-mesh2['elempartpts'].flatten('F')));
    diff[12] = mymax(abs(mesh1['cgelcon'].flatten('F')-mesh2['cgelcon'].flatten('F')));
    diff[13] = mymax(abs(mesh1['rowent2elem'].flatten('F')-mesh2['rowent2elem'].flatten('F')));
    diff[14] = mymax(abs(mesh1['cgent2dgent'].flatten('F')-mesh2['cgent2dgent'].flatten('F')));
    diff[15] = mymax(abs(mesh1['colent2elem'].flatten('F')-mesh2['colent2elem'].flatten('F')));
    diff[16] = mymax(abs(mesh1['rowe2f1'].flatten('F')-mesh2['rowe2f1'].flatten('F')));
    diff[17] = mymax(abs(mesh1['cole2f1'].flatten('F')-mesh2['cole2f1'].flatten('F')));
    diff[18] = mymax(abs(mesh1['ent2ind1'].flatten('F')-mesh2['ent2ind1'].flatten('F')));
    # diff[19] = mymax(abs(mesh1['rowe2f2'].flatten('F')-mesh2['rowe2f2'].flatten('F')));
    # diff[20] = mymax(abs(mesh1['cole2f2'].flatten('F')-mesh2['cole2f2'].flatten('F')));
    # diff[21] = mymax(abs(mesh1['ent2ind2'].flatten('F')-mesh2['ent2ind2'].flatten('F')));

    print(mesh1['rowe2f2'].shape)
    print(mesh2['rowe2f2'].shape)

    return diff
