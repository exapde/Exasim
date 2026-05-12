import numpy

def readapp(fileapp):

    app = {'nsize' : []};

    tm = numpy.fromfile(open(fileapp, "rb"), dtype=numpy.float64);

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

    app['vindx'] = numpy.array([]);
    app['dae_dt'] = numpy.array([]);
    app['interfacefluxmap'] = numpy.array([]);
    app['avparam'] = numpy.array([]);
    app['wmModelIDs'] = numpy.array([]);
    app['wmBoundaries'] = numpy.array([]);
    app['wmDistances'] = numpy.array([]);

    if len(app['nsize']) >= 13:
        k1 = k2;
        k2 = k1+app['nsize'][12];
        app['vindx'] = tm[k1:k2];
    if len(app['nsize']) >= 14:
        k1 = k2;
        k2 = k1+app['nsize'][13];
        app['dae_dt'] = tm[k1:k2];
    if len(app['nsize']) >= 15:
        k1 = k2;
        k2 = k1+app['nsize'][14];
        app['interfacefluxmap'] = tm[k1:k2];
    if len(app['nsize']) >= 16:
        k1 = k2;
        k2 = k1+app['nsize'][15];
        app['avparam'] = tm[k1:k2];
    if len(app['nsize']) >= 19:
        k1 = k2;
        k2 = k1+app['nsize'][16];
        app['wmModelIDs'] = tm[k1:k2];

        k1 = k2;
        k2 = k1+app['nsize'][17];
        app['wmBoundaries'] = tm[k1:k2];

        k1 = k2;
        k2 = k1+app['nsize'][18];
        app['wmDistances'] = tm[k1:k2];

    return app
