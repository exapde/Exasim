import os
import numpy

def masternodes(porder,dim,elemtype):

    d0 = os.getcwd();
    ii = d0.find("Exasim");
    fn = d0[0:ii] + "Exasim/src/Python/Preprocessing/masternodes.bin";
    tmp = numpy.fromfile(open(fn, "r"), dtype=numpy.float64);

    ndims = numpy.int_(tmp[0]);
    k1 = 1;
    k2 = k1+(ndims);
    narrays = numpy.int_(tmp[k1:k2]);
    k1 = k2; k2 = k1 + numpy.prod(narrays);
    sz1 = numpy.reshape(numpy.int_(tmp[k1:k2]), narrays, order='F');
    k1 = k2; k2 = k1 + numpy.prod(narrays);
    sz2 = numpy.reshape(numpy.int_(tmp[k1:k2]), narrays, order='F');
    sz = numpy.reshape(sz1*sz2, [numpy.prod(narrays), 1], order='F');
    lz = numpy.cumsum(numpy.insert(sz, 0, 0, axis=0));
    n = len(lz);
    lz1 = numpy.reshape(lz[0:(n-1)], narrays, order='F');
    lz2 = numpy.reshape(lz[1:n], narrays, order='F');

    e = elemtype;
    pm1 = porder-1;
    dm1 = dim-1;

    i = 0;
    m1 = k2 + lz1[i,e,pm1,dm1];
    m2 = m1 + lz2[i,e,pm1,dm1]-lz1[i,e,pm1,dm1];
    pelem = numpy.reshape(tmp[m1:m2],[sz1[i,e,pm1,dm1],sz2[i,e,pm1,dm1]],order='F');

    i = 1;
    m1 = k2 + lz1[i,e,pm1,dm1];
    m2 = m1 + lz2[i,e,pm1,dm1]-lz1[i,e,pm1,dm1];
    telem = numpy.reshape(tmp[m1:m2],[sz1[i,e,pm1,dm1],sz2[i,e,pm1,dm1]],order='F');

    i = 2;
    m1 = k2 + lz1[i,e,pm1,dm1];
    m2 = m1 + lz2[i,e,pm1,dm1]-lz1[i,e,pm1,dm1];
    pface = numpy.reshape(tmp[m1:m2],[sz1[i,e,pm1,dm1],sz2[i,e,pm1,dm1]],order='F');

    i = 3;
    m1 = k2 + lz1[i,e,pm1,dm1];
    m2 = m1 + lz2[i,e,pm1,dm1]-lz1[i,e,pm1,dm1];
    tface = numpy.reshape(tmp[m1:m2],[sz1[i,e,pm1,dm1],sz2[i,e,pm1,dm1]],order='F');

    i = 4;
    m1 = k2 + lz1[i,e,pm1,dm1];
    m2 = m1 + lz2[i,e,pm1,dm1]-lz1[i,e,pm1,dm1];
    perm  = numpy.reshape(tmp[m1:m2],[sz1[i,e,pm1,dm1],sz2[i,e,pm1,dm1]],order='F');

    telem = numpy.int_(telem);
    tface = numpy.int_(tface);
    perm = numpy.int_(perm);

    return pelem, telem, pface, tface, perm;
