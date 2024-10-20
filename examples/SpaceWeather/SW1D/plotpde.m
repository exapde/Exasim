close all

% Specify an Exasim version to run
version = "Version0.1";

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Install/setpath.m");

% generate input files and store them in datain folder
load('inputFiles');
tfixed = 5000;

pde.soltime = pde.saveSolFreq:pde.saveSolFreq:tfixed;

% get solution from output files in dataout folder
sol = fetchsolution(pde,master,dmd);

dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);    
figure(1);clf;
for i=1:size(sol,4)
    rho = sol(:,1,:,i); 
    vr = sol(:,2,:,i)./sqrt(exp(rho));
    T = sol(:,3,:,i)./sqrt(exp(rho));
    
    plot(dgnodes(:),T(:),'linewidth',2);
    %axis([591 593 -0.5 0.5])
    hold on
    plot(dgnodes(:),vr(:),'linewidth',2);
    plot(dgnodes(:),rho(:),'linewidth',2);
    pause(.1)
    hold off

end

% % % visualize the numerical solution of the PDE model using Paraview
% pde.visscalars = {"density", 1, "velocity radial", 2, "velocity tan", 3, "temperature", 4};  % list of scalar fields for visualization
% pde.visvectors = {}; % list of vector fields for visualization
% xdg = vis(sol,pde,mesh); % visualize the numerical solution
disp("Done!");
