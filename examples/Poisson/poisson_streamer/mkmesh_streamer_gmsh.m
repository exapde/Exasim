function mesh = mkmesh_streamer_gmsh(porder, fname)

    [p,t] = gmshcall(fname, 2, 0);

    % Normalization beacuse gmsh import is weird for some reason
    xmax = max(p(:,1));
    p = p/xmax * 125;
    
    % Below is nondimensional
    boundaryexpr = {'all(p(:,2)<min(p0(:,2))+1e-6)',...
               'all(p(:,1)>max(p0(:,1))-1e-6)', ...
               'all(p(:,2)>max(p0(:,2))-1e-6)',...
               'all(p(:,1)<min(p0(:,1))+1e-6)'};
    
    % Boundaries
    % 1 Bottom electrode
    % 2 Right farfield
    % 3 Top electrode
    % 4 Symmetry
             
    elemtype = 0;
    nodetype = 1;
    mesh = mkmesh(p',t', porder, boundaryexpr, elemtype, nodetype);