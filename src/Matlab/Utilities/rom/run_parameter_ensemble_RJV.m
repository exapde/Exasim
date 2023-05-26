function [JPhi_train] = run_parameter_ensemble_RJV(parameters, snapshots_mat, Phi, pde, master, mesh, dmd)
%Given a collection of parameters, run FOMs
% For now: parameters should be a struct
npe = master.npe;
ncu = pde.ncu;
ne = pde.ne;
u_shape = [npe ncu ne];
u_flat = [npe*ncu*ne 1];
fileapp = "datain/app.bin";
n_param = length(parameters);
p_rb = size(Phi,2);
% snapshots = cell(n_param,1);
JPhi_train = zeros(size(Phi,1), p_rb, n_param);
pde.runmode = 20;
for i = 1:n_param
    disp("Running case: " + string(i));
    param = parameters{i};
    pde.physicsparam = param;
    writeapp(pde,fileapp,'native');
    mesh.udg = reshape(snapshots_mat(:,i), u_shape);
    for ii = 1:p_rb
        % Evaluate JPhi by calling Exasim's Jv against each column of Phi
        mesh.dudg = reshape(Phi(:,ii), u_shape);
        writesol(pde,mesh,master,dmd)
        runcode(pde);
        Jv = reshape(getsolution('dataout/out_Jv_test',dmd,npe),u_flat);
        JPhi_train(:,ii,i) = Jv;
    end
%     writeapp(pde,fileapp,'native');
%     runcode(pde);
%     sol = fetchsolution(pde,master,dmd);
%     snapshots{i} = sol;
end
end