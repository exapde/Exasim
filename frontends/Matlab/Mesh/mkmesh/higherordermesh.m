function msh = higherordermesh(msh,porder)

if porder<=msh.porder
    fprintf('WARNING : higherordermesh.m is actually used to create a lower or equal order mesh.\n')
end

% if porder>msh.porder
    % higher order master nodes
    [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface] = mkmasternodes(porder,msh.nd,msh.elemtype,msh.nodetype);
    
    % shape functions and derivatives on the master volume element at higher order nodes  
    shapvl = mkshape(msh.porder,msh.plocal,plocal,msh.elemtype);
    shapvl = shapvl(:,:,1)';
    npv = size(shapvl,1);
    
    % higher order DG nodes
    dgnodes = zeros(npv, size(msh.dgnodes,2), msh.ne);
    for i = 1:msh.ne
        dgnodes(:,:,i) = shapvl*msh.dgnodes(:,:,i);
    end
    
    msh.porder = porder;
    msh.dgnodes = dgnodes;
    msh.plocal = plocal;
    msh.tlocal = tlocal;
    msh.plocfc = plocfc;
    msh.tlocfc = tlocfc;
    msh.permnode = permnode; % positions of nodes at the vertices
    msh.permedge = permedge; % positions of nodes on the edges
    msh.permface = permface; % positions of nodes on the faces
    dim = msh.nd;
    % positions of nodes on the faces
    if dim==1
        msh.perm = permnode;
    elseif dim==2
        msh.perm = permedge;
    elseif dim==3
        msh.perm = permface;
    end
    msh.permgeom = msh.perm;
% end




