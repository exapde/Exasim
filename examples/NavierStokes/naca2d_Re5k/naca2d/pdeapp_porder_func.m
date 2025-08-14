function [UDG1, mesh] = pdeapp_porder_func(pde, mesh, porder, UDG, project_flag)
    % % Add Exasim to Matlab search path
    % cdir = pwd(); ii = strfind(cdir, "Exasim"); ii = ii(end);
    % run(cdir(1:(ii+5)) + "/install/setpath.m");
    if nargin < 5
        project_flag = 0;
    end

    pde.porder = porder;
    pde.pgauss = 2*porder;
    
    % naca mesh
    % mesh = mkmesh_naca0012(porder,1,2);
    % make a higher order mesh

    % mesh size
    [~,mesh,master,dmd] = preprocessing(pde,mesh);
    mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,porder);    
    mesh.porder = porder;
    mesh.vdg = [];

    pde.pgauss = 2*(pde.porder);
    pde.nd = 2;
    pde.elemtype = 1;
    master = Master(pde);
    [~, ~, jac] = volgeom(master.shapent,permute(mesh.dgnodes,[1 3 2]));
    hsz = reshape(sqrt(jac),[],1,size(mesh.dgnodes,3));
    [~,cgelcon,rowent2elem,colent2elem,~] = mkcgent2dgent(mesh.dgnodes,1e-8);
    hh = dg2cg2(max(hsz,0e-5), cgelcon, colent2elem, rowent2elem);
    hh = dg2cg2(hh, cgelcon, colent2elem, rowent2elem);
    mesh.vdg(:,1,:) = hh;

    % if nargin == 4
    %     UDG = dgprojection(master, mesh, Uout1, pde.porder-1);
    %     % mesh.udg = UDG;
    % end
    if nargin > 3
        % UDG = dgprojection(master, mesh, Uout1, pde.porder-1);
        if project_flag
            mesh.udg = dgprojection(master, mesh, UDG, porder-1);
        else
            mesh.udg = UDG;
        end
    end


    % call exasim to generate and run C++ code to solve the PDE model
    [pde,mesh,master,dmd] = preprocessing(pde,mesh);
    if pde.gencode==1
        %gencode(pde);
        kkgencode(pde);
        compilerstr = cmakecompile(pde); % use cmake to compile C++ source codes 
    end
    runcode(pde, 1);
    %% plot solution
    % Output is saved as a steady solution (out_np0.bin) since we do not know how many timesteps are taken. 
    % So, we use getsolution; fetchsolution will grab the last timestep
    sol = getsolution(pde.buildpath + "/dataout/out", dmd, master.npe);

    UDG1 = sol;
end