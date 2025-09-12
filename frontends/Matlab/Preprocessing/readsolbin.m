function sol = readsolbin(filename)
    endian = 'native';
    fileID = fopen(filename, 'r');

    % Read number of entries in nsize
    numsize = fread(fileID, 1, 'double', endian);

    % Read size arrays
    nsize = fread(fileID, numsize, 'double', endian);
    ndims = fread(fileID, nsize(1), 'double', endian);

    % Assign basic dimensions
    sol.nsize = nsize;
    sol.ndims = ndims;

    % Read and reshape each field
    sol.xdg  = fread(fileID, nsize(2), 'double', endian);
    sol.udg  = fread(fileID, nsize(3), 'double', endian);
    sol.vdg  = fread(fileID, nsize(4), 'double', endian);
    sol.wdg  = fread(fileID, nsize(5), 'double', endian);
    sol.uhat = fread(fileID, nsize(6), 'double', endian);

    fclose(fileID);
end
