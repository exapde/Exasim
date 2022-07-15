% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";       % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel"; % name of a file defining the PDE model

% Set discretization parameters, physical parameters, and solver parameters
pde.physicsparam = [1 0.0]; % unit thermal conductivity and zero boundary data
pde.tau = 1.0;              % DG stabilization parameter
pde.NLtol = 1e-10;

% Choose computing platform and set number of processors
pde.platform = "gpu";      % choose this option if you want to run the C++ code on Nvidia GPUs
pde.mpiprocs = 1;           % number of MPI processors

ue = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2)).*sin(pi*x(:,3));
qe = @(x) pi*[cos(pi*x(:,1)).*sin(pi*x(:,2)).*sin(pi*x(:,3)), cos(pi*x(:,2)).*sin(pi*x(:,3)).*sin(pi*x(:,1)), cos(pi*x(:,3)).*sin(pi*x(:,1)).*sin(pi*x(:,2))];

nOrders = 4;
nMeshes = 5;

ErrorU = zeros(nOrders,nMeshes);
ErrorQ = zeros(nOrders,nMeshes);

ErrorU2 = zeros(nOrders,nMeshes);
ErrorQ2 = zeros(nOrders,nMeshes);


for iOrder = 1:nOrders
    for jMesh = 1:nMeshes
        mesh = initializemesh("src");
        
        pde.porder = iOrder;             % polynomial degree
        N = 2^(jMesh);
        [mesh.p,mesh.t] = cubemesh(N,N,N,0);
        % expressions for domain boundaries
        mesh.boundaryexpr = {@(p) abs(p(2,:))<1e-8, @(p) abs(p(1,:)-1)<1e-8, @(p) abs(p(2,:)-1)<1e-8, @(p) abs(p(1,:))<1e-8, @(p) abs(p(3,:))<1e-8, @(p) abs(p(3,:)-1)<1e-8};
        mesh.boundarycondition = [1;1;1;1;1;1]; % Set boundary condition for each boundary

        % call exasim to generate and run C++ code to solve the PDE model
        pde = setcompilers(pde);       
        [pde,mesh,master,dmd] = preprocessing(pde,mesh);
        gencode(pde);
        % compile source codes to build an executable file and store it in app folder
        compilerstr = compilecode(pde);

        % run executable file to compute solution and store it in dataout folder
        runstr = runcode(pde);

        % get solution from output files in dataout folder
        sol = fetchsolution(pde,master,dmd, 'dataout');
        % get residual norms from output files in dataout folder
        if pde.saveResNorm
            fileID = fopen('dataout/out_residualnorms0.bin','r'); res = fread(fileID,'double'); fclose(fileID);
            res = reshape(res,4,[])';
        end    
        
        dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);
        x = dgnodes(:,1,:); y = dgnodes(:,2,:); z = dgnodes(:,3,:);
        uexact = sin(pi*x).*sin(pi*y).*sin(pi*z);           % exact solution
        uh = sol(:,1,:);
        qh = sol(:,2:4,:);
        fprintf('Maximum absolute error: %g\n',max(abs(uh(:)-uexact(:))));
        
        mesh.dgnodes = dgnodes;
        [erru,rerru] = computeerror(uh,mesh,master,ue);
        [errq,rerrq] = computeerror(-qh,mesh,master,qe);
        ErrorU(iOrder,jMesh) = erru;
        ErrorU2(iOrder,jMesh) = rerru;
        errq = norm(errq);
        rerrq = norm(rerrq);
        ErrorQ(iOrder,jMesh) = errq;
        ErrorQ2(iOrder,jMesh) = rerrq;
        
        fprintf('L2 error u: %g\n',erru);
        fprintf('L2 relative error u: %g\n',rerru);
        fprintf('L2 error q: %g\n',errq);
        fprintf('L2 relative error q: %g\n',rerrq);
        
        rmdir dataout s
        rmdir datain s
        rmdir app s
    end
end

