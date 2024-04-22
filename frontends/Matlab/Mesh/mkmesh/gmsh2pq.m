function [p,q] = gmsh2pq(fname)

    fid = fopen(fname,'r');
    if fid == -1
        error('Can''t open file');
    end

    readuntil(fid, '$Nodes');
    np = fscanf(fid, '%d', 1);
    p = zeros(np,2);
    for ii = 1:np
        foo = fscanf(fid, '%d', 1);
        p(ii,:) = fscanf(fid, '%f', 2)';
        fgetl(fid);
    end
    
    readuntil(fid, '$Elements');
    nq = fscanf(fid, '%d', 1);
    q = zeros(nq,4);
    for ii = 1:nq
        foo = fscanf(fid, '%d', 1);
        eltype = fscanf(fid, '%d', 1);
        if eltype == 3
            foo = fscanf(fid, '%d', 3);
            q(ii,:) = fscanf(fid, '%d', 4)';
        end
        fgetl(fid);
    end
    
    q = q(any(q~=0,2),:);
    
    fclose(fid);
end

function readuntil(fid,str)
    
    while ~feof(fid)
        fline = fgetl(fid);
        if ~isempty(fline) & ~isempty(strmatch(str, fline))
            break;
        end
    end
end
