function master = readmasterbin(filename)
    endian = 'native';
    fileID = fopen(filename, 'r');

    % Read number of entries in nsize
    numsize = fread(fileID, 1, 'double', endian);

    % Read size arrays
    nsize = fread(fileID, numsize, 'double', endian);
    ndims = fread(fileID, nsize(1), 'double', endian);
    master.nsize = nsize;
    master.ndims = ndims;
    
    % Read arrays in order
    master.shapegt    = fread(fileID, nsize(2),  'double', endian);
    master.shapegw    = fread(fileID, nsize(3),  'double', endian);
    master.shapfgt    = fread(fileID, nsize(4),  'double', endian);
    master.shapfgw    = fread(fileID, nsize(5),  'double', endian);
    master.shapent    = fread(fileID, nsize(6),  'double', endian);
    master.shapen     = fread(fileID, nsize(7),  'double', endian);
    master.shapfnt    = fread(fileID, nsize(8),  'double', endian);
    master.shapfn     = fread(fileID, nsize(9),  'double', endian);
    master.xpe        = fread(fileID, nsize(10), 'double', endian);
    master.gpe        = fread(fileID, nsize(11), 'double', endian);
    master.gwe        = fread(fileID, nsize(12), 'double', endian);
    master.xpf        = fread(fileID, nsize(13), 'double', endian);
    master.gpf        = fread(fileID, nsize(14), 'double', endian);
    master.gwf        = fread(fileID, nsize(15), 'double', endian);
    master.shap1dgt   = fread(fileID, nsize(16), 'double', endian);
    master.shap1dgw   = fread(fileID, nsize(17), 'double', endian);
    master.shap1dnt   = fread(fileID, nsize(18), 'double', endian);
    master.shap1dn    = fread(fileID, nsize(19), 'double', endian);
    master.xp1d       = fread(fileID, nsize(20), 'double', endian);
    master.gp1d       = fread(fileID, nsize(21), 'double', endian);
    master.gw1d       = fread(fileID, nsize(22), 'double', endian);

    fclose(fileID);
end
