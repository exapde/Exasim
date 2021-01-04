from vtuwrite import vtuwrite
from numpy import matmul, reshape
import sys

def pvdwrite(filename, cgnodes, cgelcon, cgcells, celltype, scalars, vectors, fields, dt, visshape):

    if sys.byteorder=='little':
        byte_order = "LittleEndian";
    elif sys.byteorder=='big':
        byte_order = "BigEndian";
    else:
        sys.exit("Endian is not valid");

    pvdfile = filename + ".pvd";
    fid = open(pvdfile, "w");
    mystr = "<?xml version=\"1.0\"?>\n";
    mystr = mystr + "<VTKFile type=\"Collection\" version=\"0.1\"\n";
    mystr = mystr + "         byte_order=\"" + byte_order + "\"\n";
    mystr = mystr + "         compressor=\"vtkZLibDataCompressor\">\n";
    mystr = mystr + "  <Collection>\n";

    if fields.ndim==3:
        nt==1;
    else:
        nt = fields.shape[3];

    for i in range(0,nt):
        tm = matmul(visshape,reshape(fields[:,:,:,i],(visshape.shape[1], fields.shape[1]*fields.shape[2]), 'F'));
        tm = reshape(tm,(visshape.shape[0], fields.shape[1], fields.shape[2]),'F');

        vtufile = filename + str(i+1);
        vtuwrite(vtufile, cgnodes, cgelcon, cgcells, celltype, scalars, vectors, tm);

        ind = vtufile.find("/");
        outfile = vtufile[(ind+1):] + ".vtu";

        mystr = mystr + "    <DataSet timestep=\"" + str((i+1)*dt) + "\" group=\"\" part=\"0\"\n";
        mystr = mystr + "             file=\"" + outfile + "\"/>\n"

    fid.write(mystr);

    fid.write("  </Collection>\n");
    fid.write("</VTKFile>\n");
    fid.close();

    return 0;
