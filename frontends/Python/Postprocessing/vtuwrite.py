import sys
from numpy import *

def vtuwrite(filename, cgnodes, cgelcon, cgcells, celltype, scalars, vectors, fields):

    cgelcon = cgelcon.flatten('F');
    filename = filename + ".vtu";

    # get dimensions
    npoints,nd = shape(cgnodes);
    ncells,nve = shape(cgcells);

    # Preparing output for CG
    if nd==2:
        cgnodes = hstack([cgnodes,zeros((npoints,1))]);
    outputCG = zeros((npoints,3));

    # Get endianness
    if sys.byteorder=='little':
        byte_order = "LittleEndian";
    elif sys.byteorder=='big':
        byte_order = "BigEndian";
    else:
        sys.exit("Endian is not valid");

    # binary format
    formattype = "appended";

    # float and integer byte size
    fbytesize = 4; # Float32 format
    ibytesize = 4; # Int32   format
    obytesize = 8; # Int64   format

    nsc = len(scalars);
    if nsc>1:
        if nsc==2:
            sidx = [0];
            sidy = [1];
        else:
            sidx = arange(0,nsc,2);
            sidy = arange(1,nsc+1,2);
    else:
        sidx = [];
        sidy = [];

    nvt = len(vectors);
    if nvt>1:
        if nvt==2:
            vidx = [0];
            vidy = [1];
        else:
            vidx = arange(0,nvt,2);
            vidy = arange(1,nvt+1,2);
    else:
        vidx = [];
        vidy = [];

    # Open VTK output file
    fid = open(filename, "w");
    # VTK DataFile Version
    mystr = "<?xml version=\"1.0\"?>\n";
    mystr = mystr + "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"" + byte_order + "\" header_type=\"UInt64\">\n";
    mystr = mystr + "  <UnstructuredGrid>\n";
    mystr = mystr + "    <Piece NumberOfPoints=\"" + str(npoints) + "\"" +  " NumberOfCells=\"" + str(ncells) + "\">\n";

    ###########################################################################
    #                      1 - WRITE VTU METADATA                             #
    #    (Write all dataset description for the VTU binary-appended format    #
    ###########################################################################
    offset = 0;
    if (len(sidx)>0) or (len(vidx)>0):
        mystr = mystr + "      <PointData Scalars=\"scalars\">\n";

    for ii in range(0,len(sidx)):
        title = scalars[sidx[ii]];
        mystr = mystr + "        <DataArray type=\"Float32\" Name=\"" + title + "\" Format=\"" + formattype + "\" offset=\"" + str(offset) + "\">\n";
        mystr = mystr + "        </DataArray>\n";
        offset = offset + npoints*fbytesize + obytesize;

    for ii in range(0,len(vidx)):
        title = vectors[vidx[ii]];
        mystr = mystr + "        <DataArray type=\"Float32\" Name=\"" + title + "\" NumberOfComponents=\"3" + "\" Format=\"" + formattype + "\" offset=\"" + str(offset) + "\">\n";
        mystr = mystr + "        </DataArray>\n";
        offset = offset + 3*npoints*fbytesize + obytesize;

    if (len(sidx)>0) or (len(vidx)>0):
        mystr = mystr + "      </PointData>\n";

    mystr = mystr + "      <Points>\n";
    mystr = mystr + "        <DataArray type=\"Float32\" Name=\"" + "points" + "\" NumberOfComponents=\"3" + "\" Format=\"" + formattype + "\" offset=\"" + str(offset) + "\">\n";
    mystr = mystr + "        </DataArray>\n";
    mystr = mystr + "      </Points>\n";
    offset = offset + 3*npoints*fbytesize + obytesize;
    mystr = mystr + "      <Cells>\n";
    mystr = mystr + "        <DataArray type=\"Int32\" Name=\"" + "connectivity" + "\" Format=\"" + formattype + "\" offset=\"" + str(offset) + "\">\n";
    mystr = mystr + "        </DataArray>\n";
    offset = offset + ncells*nve*ibytesize + obytesize;
    mystr = mystr + "        <DataArray type=\"Int32\" Name=\"" + "offsets" + "\" Format=\"" + formattype + "\" offset=\"" + str(offset) + "\">\n";
    mystr = mystr + "        </DataArray>\n";
    offset = offset + ncells*ibytesize + obytesize;
    mystr = mystr + "        <DataArray type=\"UInt8\" Name=\"" + "types" + "\" Format=\"" + formattype + "\" offset=\"" + str(offset) + "\">\n";
    mystr = mystr + "        </DataArray>\n";
    mystr = mystr + "      </Cells>\n";
    mystr = mystr + "    </Piece>\n";
    mystr = mystr + "  </UnstructuredGrid>\n";
    mystr = mystr + "  <AppendedData encoding=\"raw\">\n"
    mystr = mystr + "   _";

    fid.write(mystr);

    bsize = zeros(1).astype(int64);

    ###########################################################################
    #                   2 - WRITE UNSTRUCTURED DATA                           #
    #        (Vectors and Scalars datasets all projected on CG mesh           #
    ###########################################################################
    # This part writes the dataset attributes for the mesh nodes or elements.
    # It is customary to write the dataset before the mesh description.

    #  Write all the scalar dataset first
    for ii in range(0,len(sidx)):
        output = fields[:,scalars[sidy[ii]],:];#varargin{ sidx(ii) + 2};
        outputCG[cgelcon,0] = output.flatten('F');
        bsize[0] = npoints*fbytesize;
        bsize.astype('int64').tofile(fid);
        tmp = outputCG[:,0];
        tmp.astype('float32').tofile(fid);

    #  Write all the vector datasets then
    for ii in range(0,len(vidx)):
        tmp = fields[:,vectors[vidy[ii]][0],:];
        output = zeros((size(tmp),3));
        output[:,0] = tmp.flatten('F');
        tmp = fields[:,vectors[vidy[ii]][1],:];
        output[:,1] = tmp.flatten('F');
        if nd==3:
            tmp = fields[:,vectors[vidy[ii]][2],:];
            output[:,2] = tmp.flatten('F');

        # Convert to CG Field
        outputCG[cgelcon,0] = output[:,0];
        outputCG[cgelcon,1] = output[:,1];
        outputCG[cgelcon,2] = output[:,2];
        # For a more accurate conversion, use the slower
        bsize[0] = 3*npoints*fbytesize;
        bsize.astype('int64').tofile(fid);
        #tmp = outputCG.transpose().flatten('F');
        #tmp.astype('float32').tofile(fid);
        outputCG.astype('float32').tofile(fid);

    ###########################################################################
    #                3 - WRITE UNSTRUCTURED_GRID MESH                         #
    #       (Mesh nodes coordinates connectivities, offset ant type)          #
    ###########################################################################
    # Write 3D coordinates of mesh points
    bsize[0] = 3*npoints*fbytesize;
    bsize.astype('int64').tofile(fid);
    # output = cgnodes.transpose().flatten('F');
    # output.astype('float32').tofile(fid);
    cgnodes.astype('float32').tofile(fid);

    # Write (sub)-cells connectivities
    cgcells = cgcells.astype(int32);
    bsize[0] = nve*ncells*ibytesize;
    bsize.astype('int64').tofile(fid);
    # output = cgcells.transpose().flatten('F');
    # output.astype('int32').tofile(fid);
    cgcells.astype('int32').tofile(fid);

    # Write (sub)-cells offsets
    bsize[0] = ncells*ibytesize;
    bsize.astype('int64').tofile(fid);
    output = arange(nve, nve*(ncells+1), nve);
    output.astype('int32').tofile(fid);

    # Write (sub)-cells type
    bsize[0] = ncells;
    bsize.astype('int64').tofile(fid);
    output = celltype*(ones(ncells).astype('uint8'));
    output.astype('uint8').tofile(fid);

    fid.write("  </AppendedData>\n");
    fid.write("</VTKFile>\n");
    fid.close();

    return 0;
