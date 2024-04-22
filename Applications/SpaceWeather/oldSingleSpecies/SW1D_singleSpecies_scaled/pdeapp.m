% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();
pde.saveResNorm = 1;

% Choose computing platform and set number of processors
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

[pde,mesh] = pdeparams(pde,mesh);

% search compilers and set options
pde = setcompilers(pde);       

% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

save('inputFiles','pde','mesh','master','dmd');

% generate source codes and store them in app folder
gencode(pde);

% compile source codes to build an executable file and store it in app folder
compilerstr = compilecode(pde);

% % run executable file to compute solution and store it in dataout folder
runstr = runcode(pde);

%get solution from output files in dataout folder
sol = fetchsolution(pde,master,dmd);
res = fetchresidual(pde);

dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);    
figure(1);clf;
for i=1:size(sol,4)
    rho = sol(:,1,:,i); 
    vr = sol(:,2,:,i)./sqrt(exp(rho));
    T = sol(:,3,:,i)./sqrt(exp(rho));
    
    plot(dgnodes(:),T(:),'linewidth',2);
    axis([min(dgnodes(:)) max(dgnodes(:)) -18 10])
    hold on
    grid on
    grid minor
    plot(dgnodes(:),vr(:),'linewidth',2);
    plot(dgnodes(:),rho(:),'linewidth',2);
    
    text(min(dgnodes(:)) + 10,-11,sprintf('t = %d',i))

    
    pause(.1)
    hold off
end

disp("Done!");
