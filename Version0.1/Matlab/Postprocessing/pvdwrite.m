function pvdwrite(filename, cgnodes, cgelcon, cgcells, celltype, scalars, vectors, fields, dt)

[~,~,endian] = computer;
if strcmp(endian,'L')
    byte_order = 'LittleEndian';
else
    byte_order = 'BigEndian';
end

% Write Header for Paraview PVD File
pvdfile = filename + ".pvd";
fid = fopen(pvdfile, 'w');
fprintf(fid,'<?xml version="1.0"?>\n');
fprintf(fid,'<VTKFile type="Collection" version="0.1"\n');
fprintf(fid,'         byte_order="%s"\n',byte_order);
fprintf(fid,'         compressor="vtkZLibDataCompressor">\n');
fprintf(fid,'  <Collection>\n');

    
nt = size(fields,4);
for i = 1:nt
    vtufile = filename + num2str(i);
    
    vtuwrite(vtufile, cgnodes, cgelcon, cgcells, celltype, scalars, vectors, fields(:,:,:,i));
        
    ind = strfind(vtufile, "/");
    ch = char(vtufile);
    outfile = string(ch((ind(end)+1):end));
    
    % Write pointer to .vtu file in .pvd file
    fprintf(fid,'    <DataSet timestep="%e" group="" part="0"\n',i*dt);
    fprintf(fid,'             file="%s"/>\n', outfile + ".vtu");    
end

% Finalize and close .pvd file
fprintf(fid,'  </Collection>\n');
fprintf(fid,'</VTKFile>\n');
fclose(fid);

