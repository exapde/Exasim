import numpy

def readmeshstruct(filemesh):

    mesh = {'nsize' : []};

    tm = numpy.fromfile(open(filemesh, "r"), dtype=numpy.float64);
    tm = numpy.int_(tm);

    sz = numpy.int_(tm[0]);
    k1 = 1;
    k2 = k1+(sz);
    mesh['nsize'] = numpy.int_(tm[k1:k2]);

    k1 = k2;
    k2 = k1+mesh['nsize'][0];
    mesh['ndims'] = numpy.int_(tm[k1:k2]);

    k1 = k2;
    k2 = k1+mesh['nsize'][1];
    mesh['facecon'] = numpy.int_(tm[k1:k2]);

    k1 = k2;
    k2 = k1+mesh['nsize'][2];
    mesh['eblks'] = numpy.int_(tm[k1:k2]);

    k1 = k2;
    k2 = k1+mesh['nsize'][3];
    mesh['fblks'] = numpy.int_(tm[k1:k2]);

    k1 = k2;
    k2 = k1+mesh['nsize'][4];
    mesh['nbsd'] = numpy.int_(tm[k1:k2]);

    k1 = k2;
    k2 = k1+mesh['nsize'][5];
    mesh['elemsend'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][6];
    mesh['elemrecv'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][7];
    mesh['elemsendpts'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][8];
    mesh['elemrecvpts'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][9];
    mesh['elempart'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][10];
    mesh['elempartpts'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][11];
    mesh['cgelcon'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][12];
    mesh['rowent2elem'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][13];
    mesh['cgent2dgent'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][14];
    mesh['colent2elem'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][15];
    mesh['rowe2f1'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][16];
    mesh['cole2f1'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][17];
    mesh['ent2ind1'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][18];
    mesh['rowe2f2'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][19];
    mesh['cole2f2'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+mesh['nsize'][20];
    mesh['ent2ind2'] = tm[k1:k2];

    return mesh
