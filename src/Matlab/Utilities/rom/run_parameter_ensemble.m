function [snapshots] = run_parameter_ensemble(parameters, pde, master, dmd)
%Given a collection of parameters, run FOMs
% For now: parameters should be a struct
fileapp = "datain/app.bin";
n_param = length(parameters);
snapshots = cell(n_param,1);
pde.runmode = 0;
for i = 1:n_param
    disp("Running case: " + string(i));
    param = parameters{i};
    pde.physicsparam = param;
    writeapp(pde,fileapp,'native');
    runcode(pde);
    sol = fetchsolution(pde,master,dmd);
    snapshots{i} = sol;
end
end