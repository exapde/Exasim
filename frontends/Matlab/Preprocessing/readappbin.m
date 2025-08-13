function app = readappbin(filename)
    endian = 'native';
    fileID = fopen(filename, 'r');
    
    % Read number of entries in nsize
    numsize = fread(fileID, 1, 'double', endian);
    
    % Read size arrays
    nsize = fread(fileID, numsize, 'double', endian);
    ndims = fread(fileID, nsize(1), 'double', endian);
    
    % Assign basic dimensions to app struct
    app.nsize = nsize;
    app.ndims = ndims;
    
    % Read standard fields
    app.flag           = fread(fileID, nsize(2),  'double', endian);
    app.problem        = fread(fileID, nsize(3),  'double', endian);
    app.externalparam  = fread(fileID, nsize(4),  'double', endian);
    app.dt             = fread(fileID, nsize(5),  'double', endian);
    app.factor         = fread(fileID, nsize(6),  'double', endian);
    app.physicsparam   = fread(fileID, nsize(7),  'double', endian);
    app.solversparam   = fread(fileID, nsize(8),  'double', endian);
    app.tau            = fread(fileID, nsize(9),  'double', endian);
    app.stgdata        = fread(fileID, nsize(10), 'double', endian);
    app.stgparam       = fread(fileID, nsize(11), 'double', endian);
    app.stgib          = fread(fileID, nsize(12), 'double', endian);    
    app.vindx = fread(fileID, nsize(13), 'double', endian);
    app.dae_dt = fread(fileID, nsize(14), 'double', endian);
    app.interfacefluxmap = fread(fileID, nsize(15), 'double', endian);
    app.avparam = fread(fileID, nsize(16), 'double', endian);
    
    fclose(fileID);
end
