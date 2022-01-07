import numpy

def readapp(fileapp):

    app = {'nsize' : []};

    tm = numpy.fromfile(open(fileapp, "r"), dtype=numpy.float64);

    sz = numpy.int_(tm[0]);
    k1 = 1;
    k2 = k1+(sz);
    app['nsize'] = numpy.int_(tm[k1:k2]);

    k1 = k2;
    k2 = k1+app['nsize'][0];
    app['ndims'] = numpy.int_(tm[k1:k2]);

    k1 = k2;
    k2 = k1+app['nsize'][1];
    app['flag'] = numpy.int_(tm[k1:k2]);

    k1 = k2;
    k2 = k1+app['nsize'][2];
    app['problem'] = numpy.int_(tm[k1:k2]);

    k1 = k2;
    k2 = k1+app['nsize'][3];
    app['uinf'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+app['nsize'][4];
    app['dt'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+app['nsize'][5];
    app['factor'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+app['nsize'][6];
    app['physicsparam'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+app['nsize'][7];
    app['solversparam'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+app['nsize'][8];
    app['tau'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+app['nsize'][9];
    app['stgdata'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+app['nsize'][10];
    app['stgparam'] = tm[k1:k2];

    k1 = k2;
    k2 = k1+app['nsize'][11];
    app['stgib'] = tm[k1:k2];

    return app
