function [p, t, h, ob, lb] = readfort14(filename)
%READFORT14 Read an ADCIRC fort.14 mesh file.
%
%   [p, t, h, ob, lb] = readfort14(filename)
%
% Inputs:
%   filename : string
%       Path to ADCIRC fort.14 file
%
% Outputs:
%   p  : np x 2 double
%       Node coordinates [x, y]
%
%   t  : ne x 3 int
%       Triangle connectivity (1-based indices as stored in fort.14)
%
%   h  : np x 1 double
%       Bathymetric depth at nodes
%
%   ob : cell array
%       Open boundaries. Each ob{i} is a column vector of node indices.
%
%   lb : cell array
%       Land boundaries. Each lb{i} is a struct with fields:
%           .nodes  : column vector of node indices
%           .ibtype : boundary type
%
% Notes:
%   1. This function assumes a standard 2D ADCIRC fort.14 format.
%   2. Depth values are returned exactly as stored in fort.14.
%      If you want positive bathymetric depth h = -z, adjust as needed.
%   3. Connectivity t is returned with MATLAB 1-based indexing.

    fid = fopen(filename, 'r');
    if fid < 0
        error('Could not open file: %s', filename);
    end

    cleaner = onCleanup(@() fclose(fid));

    %% --------------------------------------------------------------------
    %  Header
    %% --------------------------------------------------------------------
    gridname = fgetl(fid); %#ok<NASGU>
    if ~ischar(gridname)
        error('Invalid fort.14 file: missing title line.');
    end

    vals = sscanf(fgetl(fid), '%d %d');
    if numel(vals) < 2
        error('Invalid fort.14 file: missing NE and NP.');
    end
    ne = vals(1);
    np = vals(2);

    %% --------------------------------------------------------------------
    %  Nodes
    %  Format:
    %    node_id   x   y   depth
    %% --------------------------------------------------------------------
    p = zeros(np, 2);
    h = zeros(np, 1);

    for i = 1:np
        line = fgetl(fid);
        if ~ischar(line)
            error('Unexpected end of file while reading node data.');
        end

        vals = sscanf(line, '%f');
        if numel(vals) < 4
            error('Invalid node record at node %d.', i);
        end

        node_id = vals(1);
        if node_id < 1 || node_id > np
            error('Node id out of range at node record %d.', i);
        end

        p(node_id, 1) = vals(2);
        p(node_id, 2) = vals(3);
        h(node_id)    = vals(4);
    end

    %% --------------------------------------------------------------------
    %  Elements
    %  Format:
    %    elem_id   3   n1   n2   n3
    %% --------------------------------------------------------------------
    t = zeros(ne, 3);

    for i = 1:ne
        line = fgetl(fid);
        if ~ischar(line)
            error('Unexpected end of file while reading element data.');
        end

        vals = sscanf(line, '%d');
        if numel(vals) < 5
            error('Invalid element record at element %d.', i);
        end

        elem_id = vals(1);
        nverts  = vals(2);

        if nverts ~= 3
            error('Element %d has %d vertices. Only triangles are supported.', elem_id, nverts);
        end

        if elem_id < 1 || elem_id > ne
            error('Element id out of range at element record %d.', i);
        end

        t(elem_id, :) = vals(3:5).';
    end

    %% --------------------------------------------------------------------
    %  Open boundaries
    %
    %  Format:
    %    NOPE
    %    NETA
    %    For each open boundary:
    %        NVDLL
    %        node_1
    %        ...
    %        node_NVDLL
    %% --------------------------------------------------------------------
    line = fgetl(fid);
    if ~ischar(line)
        error('Unexpected end of file while reading number of open boundaries.');
    end
    nope = sscanf(line, '%d', 1);

    line = fgetl(fid);
    if ~ischar(line)
        error('Unexpected end of file while reading total open boundary nodes.');
    end
    neta = sscanf(line, '%d', 1); %#ok<NASGU>

    ob = cell(nope, 1);

    for ib = 1:nope
        line = fgetl(fid);
        if ~ischar(line)
            error('Unexpected end of file while reading open boundary header %d.', ib);
        end

        nvdll = sscanf(line, '%d', 1);
        if isempty(nvdll) || nvdll < 1
            error('Invalid open boundary header at boundary %d.', ib);
        end

        nodes = zeros(nvdll, 1);
        for k = 1:nvdll
            line = fgetl(fid);
            if ~ischar(line)
                error('Unexpected end of file while reading open boundary %d.', ib);
            end
            node = sscanf(line, '%d', 1);
            if isempty(node)
                error('Invalid open boundary node in boundary %d.', ib);
            end
            nodes(k) = node;
        end

        ob{ib} = nodes;
    end

    %% --------------------------------------------------------------------
    %  Land boundaries
    %
    %  Format:
    %    NBOU
    %    NVEL
    %    For each land boundary:
    %        NVELL IBTYPE
    %        node_1
    %        ...
    %        node_NVELL
    %
    %  Note:
    %    Some ADCIRC variants append extra data on each boundary-node line
    %    for certain IBTYPE values. This reader only extracts the first
    %    integer on each such line, which is the node id.
    %% --------------------------------------------------------------------
    line = fgetl(fid);
    if ~ischar(line)
        error('Unexpected end of file while reading number of land boundaries.');
    end
    nbou = sscanf(line, '%d', 1);

    line = fgetl(fid);
    if ~ischar(line)
        error('Unexpected end of file while reading total land boundary nodes.');
    end
    nvel = sscanf(line, '%d', 1); %#ok<NASGU>

    lb = cell(nbou, 1);

    for ib = 1:nbou
        line = fgetl(fid);
        if ~ischar(line)
            error('Unexpected end of file while reading land boundary header %d.', ib);
        end

        vals = sscanf(line, '%d');
        if numel(vals) < 2
            error('Invalid land boundary header at boundary %d.', ib);
        end

        nvell  = vals(1);
        ibtype = vals(2);

        if nvell < 1
            error('Invalid number of nodes in land boundary %d.', ib);
        end

        nodes = zeros(nvell, 1);
        for k = 1:nvell
            line = fgetl(fid);
            if ~ischar(line)
                error('Unexpected end of file while reading land boundary %d.', ib);
            end

            vals = sscanf(line, '%f');
            if isempty(vals)
                error('Invalid land boundary node in boundary %d.', ib);
            end

            nodes(k) = vals(1);
        end

        s.nodes  = nodes;
        s.ibtype = ibtype;
        lb{ib}   = s;
    end
end