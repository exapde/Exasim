function pvdwrite(filename, cgnodes, cgelcon, cgcells, celltype, scalars, vectors, fields, dt, shape)

if ENDIAN_BOM==0x04030201
    byte_order = "LittleEndian";
elseif ENDIAN_BOM==0x01020304
    byte_order = "BigEndian";
else
    error("Endian is not valid");
end

pvdfile = filename * ".pvd";
fid = open(pvdfile, "w");
mystr = "<?xml version=\"1.0\"?>\n";
mystr = mystr * "<VTKFile type=\"Collection\" version=\"0.1\"\n";
mystr = mystr * "         byte_order=\"" * byte_order * "\"\n";
mystr = mystr * "         compressor=\"vtkZLibDataCompressor\">\n";
mystr = mystr * "  <Collection>\n";

nt = size(fields,4);
for i = 1:nt
    vtufile = filename * string(i);

    ind = findlast("/", vtufile);
    outfile = vtufile[(ind[end]+1):end] * ".vtu";

    tm = shape*reshape(fields[:,:,:,i],(size(shape,2), size(fields,2)*size(fields,3)));
    tm = reshape(tm,(size(shape,1), size(fields,2), size(fields,3)));

    vtuwrite(vtufile, cgnodes, cgelcon, cgcells, celltype, scalars, vectors, tm);
    mystr = mystr * "    <DataSet timestep=\"" * string(i*dt) * "\" group=\"\" part=\"0\"\n";
    mystr = mystr * "             file=\"" * outfile * "\"/>\n"
end
write(fid,mystr);

write(fid,"  </Collection>\n");
write(fid,"</VTKFile>\n");
close(fid);

end
