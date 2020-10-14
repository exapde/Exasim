function vtuwrite(filename, cgnodes, cgelcon, cgcells, celltype, scalars, vectors, fields)

filename = filename * ".vtu";

# get dimensions
npoints,nd = size(cgnodes);
ncells,nve = size(cgcells);

if nd==2
    cgnodes = hcat(cgnodes,zeros(npoints,1));
end
outputCG = zeros(npoints,3);

if ENDIAN_BOM==0x04030201
    byte_order = "LittleEndian";
elseif ENDIAN_BOM==0x01020304
    byte_order = "BigEndian";
else
    error("Endian is not valid");
end

# binary format
formattype = "appended";

# float and integer byte size
fbytesize = 4; # Float32 format
ibytesize = 4; # Int32   format
obytesize = 8; # Int64   format

nsc = length(scalars);
if nsc>1
    if nsc==2
        sidx = [1];
        sidy = [2];
    else
        sidx = collect(1:2:(nsc-1));
        sidy = collect(2:2:nsc);
    end
else
    sidx = [];
    sidy = [];
end
nvt = length(vectors);
if nvt>1
    if nvt==2
        vidx = [1];
        vidy = [2];
    else
        vidx = collect(1:2:(nvt-1));
        vidy = collect(2:2:nvt);
    end
else
    vidx = [];
    vidy = [];
end

# Open VTK output file
fid = open(filename, "w");
# VTK DataFile Version
mystr = "<?xml version=\"1.0\"?>\n";
mystr = mystr * "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"" * byte_order * "\" header_type=\"UInt64\">\n";
mystr = mystr * "  <UnstructuredGrid>\n";
mystr = mystr * "    <Piece NumberOfPoints=\"" * string(npoints) * "\"" *  " NumberOfCells=\"" * string(ncells) * "\">\n";

offset = 0;
if (length(sidx)>0) || (length(vidx)>0)
    mystr = mystr * "      <PointData Scalars=\"scalars\">\n";
end
for ii = 1:length(sidx)
    title = scalars[sidx[ii]];
    mystr = mystr * "        <DataArray type=\"Float32\" Name=\"" * title * "\" Format=\"" * formattype * "\" offset=\"" * string(offset) * "\">\n";
    mystr = mystr * "        </DataArray>\n";
    offset = offset + npoints*fbytesize + obytesize;
end

for ii = 1:length(vidx)
    title = vectors[vidx[ii]];
    mystr = mystr * "        <DataArray type=\"Float32\" Name=\"" * title * "\" NumberOfComponents=\"3" * "\" Format=\"" * formattype * "\" offset=\"" * string(offset) * "\">\n";
    mystr = mystr * "        </DataArray>\n";
    offset = offset + 3*npoints*fbytesize + obytesize;
end
if (length(sidx)>0) || (length(vidx)>0)
    mystr = mystr * "      </PointData>\n";
end
mystr = mystr * "      <Points>\n";
mystr = mystr * "        <DataArray type=\"Float32\" Name=\"" * "points" * "\" NumberOfComponents=\"3" * "\" Format=\"" * formattype * "\" offset=\"" * string(offset) * "\">\n";
mystr = mystr * "        </DataArray>\n";
mystr = mystr * "      </Points>\n";
offset = offset + 3*npoints*fbytesize + obytesize;
mystr = mystr * "      <Cells>\n";
mystr = mystr * "        <DataArray type=\"Int32\" Name=\"" * "connectivity" * "\" Format=\"" * formattype * "\" offset=\"" * string(offset) * "\">\n";
mystr = mystr * "        </DataArray>\n";
offset = offset + ncells*nve*ibytesize + obytesize;
mystr = mystr * "        <DataArray type=\"Int32\" Name=\"" * "offsets" * "\" Format=\"" * formattype * "\" offset=\"" * string(offset) * "\">\n";
mystr = mystr * "        </DataArray>\n";
offset = offset + ncells*ibytesize + obytesize;
mystr = mystr * "        <DataArray type=\"UInt8\" Name=\"" * "types" * "\" Format=\"" * formattype * "\" offset=\"" * string(offset) * "\">\n";
mystr = mystr * "        </DataArray>\n";
mystr = mystr * "      </Cells>\n";
mystr = mystr * "    </Piece>\n";
mystr = mystr * "  </UnstructuredGrid>\n";
mystr = mystr * "  <AppendedData encoding=\"raw\">\n"
mystr = mystr * "   _";

write(fid,mystr);

#  Write all the scalar dataset first
for ii = 1:length(sidx)
    output = fields[:,scalars[sidy[ii]],:];#varargin{ sidx(ii) + 2};
    outputCG[cgelcon[:],1] = output[:];
    write(fid, Int64(npoints*fbytesize));
    write(fid, Float32.(outputCG[:,1]));
end

#  Write all the vector datasets then
for ii = 1:length(vidx)
    tmp = fields[:,vectors[vidy[ii]][1],:];
    output = zeros(length(tmp),3);
    output[:,1] = tmp[:];
    tmp = fields[:,vectors[vidy[ii]][2],:];
    output[:,2] = tmp[:];
    if nd==3
        tmp = fields[:,vectors[vidy[ii]][3],:];
        output[:,3] = tmp[:];
    end
    # Convert to CG Field
    outputCG[cgelcon[:],1] = output[:,1];
    outputCG[cgelcon[:],2] = output[:,2];
    outputCG[cgelcon[:],3] = output[:,3];
    # For a more accurate conversion, use the slower
    write(fid, Int64(3*npoints*fbytesize));
    write(fid, Float32.(outputCG'));
end

# Write 3D coordinates of mesh points
output = cgnodes';
write(fid, Int64(3*npoints*fbytesize));
write(fid, Float32.(output));

# Write (sub)-cells connectivities
output = cgcells';
write(fid, Int64(nve*ncells*ibytesize));
write(fid, Int32.(output));

# Write (sub)-cells offsets
output = collect(nve:nve:nve*ncells);
write(fid, Int64(ncells*ibytesize));
write(fid, Int32.(output));

# Write (sub)-cells type
output = celltype*ones(Int,1,ncells);
write(fid, Int64(ncells));
write(fid, UInt8.(output));

write(fid,"  </AppendedData>\n");
write(fid,"</VTKFile>\n");
close(fid);

end
