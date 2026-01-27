function [] = init_unsteady_app(app, mesh)

    app_tdep = app;

    app_tdep.GMRESrestart = 100;
    app_tdep.linearsolvertol = 1e-4; % GMRES tolerance
    app_tdep.linearsolveriter = 100;

    app_tdep.NLtol = 1e-12;              % Newton tolerance
    app_tdep.NLiter = 1;                 % Newton iterations
    app_tdep.dt = 1e-4*ones(1,1);   % time step sizes
    app_tdep.soltime = [1]; % steps at which solution are collected
    app_tdep.torder = 1;          % time-stepping order of accuracy
    app_tdep.nstage = 1;          % time-stepping number of stages
    app_tdep.visdt = app_tdep.dt(1);         % visualization timestep size
    app_tdep.saveSolFreq = 1;          % solution is saved every 10 time steps
    app_tdep.RBdim = 0;

    [app_tdep,mesh,master,dmd] = preprocessing(app_tdep,mesh);
    runstr = "!cp ./datain/app.bin ./datain/app_ptc.bin";
    eval(char(runstr));
end
