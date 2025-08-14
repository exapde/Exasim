function mesh = readmeshbin(filename)
    endian = 'native';
    fileID = fopen(filename, 'r');

    % Read number of entries in nsize
    numsize = fread(fileID, 1, 'double', endian);

    % Read nsize and ndims arrays
    nsize = fread(fileID, numsize, 'double', endian);
    ndims = fread(fileID, nsize(1), 'double', endian);

%     % Extract useful dimension info
%     nd = ndims(1);                  % dim of p
%     np = ndims(2);                  % length(elempart)
%     nfp = ndims(3);                 % length(facepartpts)
%     nt = ndims(4);                  % max index in t
%     nfe = ndims(5); mesh.nfe = nfe;
%     nbe = ndims(6); mesh.nbe = nbe;
%     neb = ndims(7); mesh.neb = neb;
%     nbf = ndims(8); mesh.nbf = nbf;
%     nfb = ndims(9); mesh.nfb = nfb;

    % Read flat arrays in the order they were written
    mesh.nsize = nsize;
    mesh.ndims = ndims;

    mesh.facecon = fread(fileID, nsize(2), 'double', endian);
    mesh.eblks = fread(fileID, nsize(3), 'double', endian);
    mesh.fblks = fread(fileID, nsize(4), 'double', endian);
    mesh.nbsd = fread(fileID, nsize(5), 'double', endian);
    mesh.elemsend = fread(fileID, nsize(6), 'double', endian);
    mesh.elemrecv = fread(fileID, nsize(7), 'double', endian);
    mesh.elemsendpts = fread(fileID, nsize(8), 'double', endian);
    mesh.elemrecvpts = fread(fileID, nsize(9), 'double', endian);
    mesh.elempart = fread(fileID, nsize(10), 'double', endian);
    mesh.elempartpts = fread(fileID, nsize(11), 'double', endian);
    mesh.cgelcon = fread(fileID, nsize(12), 'double', endian);
    mesh.rowent2elem = fread(fileID, nsize(13), 'double', endian);
    mesh.cgent2dgent = fread(fileID, nsize(14), 'double', endian);
    mesh.colent2elem = fread(fileID, nsize(15), 'double', endian);
    mesh.rowe2f1 = fread(fileID, nsize(16), 'double', endian);
    mesh.cole2f1 = fread(fileID, nsize(17), 'double', endian);
    mesh.ent2ind1 = fread(fileID, nsize(18), 'double', endian);
    mesh.rowe2f2 = fread(fileID, nsize(19), 'double', endian);
    mesh.cole2f2 = fread(fileID, nsize(20), 'double', endian);
    mesh.ent2ind2 = fread(fileID, nsize(21), 'double', endian);

    % Handle optional hybrid fields
    mesh.f2t = fread(fileID, nsize(22), 'double', endian);
    mesh.elemcon = fread(fileID, nsize(23), 'double', endian);
    mesh.perm = fread(fileID, nsize(24), 'double', endian);
    mesh.bf = fread(fileID, nsize(25), 'double', endian);
    mesh.cartgridpart = fread(fileID, nsize(26), 'double', endian);

    % Set dimensional info
%     mesh.nd = nd;
%     mesh.np = np;
%     mesh.nfp = nfp;
%     mesh.nt = nt;

    fclose(fileID);
end
