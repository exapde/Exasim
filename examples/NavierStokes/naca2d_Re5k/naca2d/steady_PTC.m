function [UDG, UH, a] = steady_PTC(app, mesh, master, UDG, UH, a)
% % setapp_dirklicationpath('FM/eulershock2');
% % % app.S0=0.02; app.lambda = 0.01*avfactor; app.kappa=1.5;
% app.S0=0.01; app.lambda = 0.01*avfactor; app.kappa=1.5;
% % 
% a = avf(mesh, master, app, UDG);
% figure(1); clf; scaplot(mesh, a,[],2); 
ns = 1;

%
app_tdep = app;
app_tdep.GMRESrestart = 40;
app_tdep.linearsolvertol = 1e-4; % GMRES tolerance
app_tdep.linearsolveriter = 80;
app_tdep.NLtol = 1e-12;              % Newton tolerance
app_tdep.NLiter = 1;                 % Newton iterations
app_tdep.dt = 1e-4*ones(1,1);   % time step sizes
app_tdep.soltime = [1]; % steps at which solution are collected
app_tdep.torder = 1;          % time-stepping order of accuracy
app_tdep.nstage = 1;          % time-stepping number of stages
app_tdep.visdt = app_tdep.dt(1);         % visualization timestep size
app_tdep.saveSolFreq = 1;          % solution is saved every 10 time steps
app_tdep.RBdim = 0;
%%
% timestep sizes

% dt = dt*0 +0.0025/2;
dt = 1e-3;
% dt = 20;
time = 0;
duh = 1;
app_tdep.read_uh = 0;

% mesh.vdg(:,1,:) = a;
for itime = 1:1000      
    revert_flag = 0;
    fprintf('Timestep :  %d,  Time :   %g\n', itime, time);
    UDG_old = UDG;
    UH_old = UH;
    %[UDG,UH] = hdg_solve_dirk(master,mesh,app_dirk,UDG,UH,[],time,dt(itime),nstage,torder);      
    %[Un,Hn,Pn] = hdg_solve_dirk(master,mesh,app_dirk,UDG,UH,PDG,time,dt,nstage,torder)
    %[Un,Hn,Pn] = hdg_solve_dirk(master,mesh,app_dirk,UDG,UH,PDG,time,dt,nstage,torder)
    mesh.udg = UDG;
    [app_tdep,mesh,master,dmd] = preprocessing(app_tdep,mesh);

    runcode(app_tdep, 1); % run C++ code
    % UDG = fetchsolution(app_tdep,master,dmd, app_tdep.buildpath + '/dataout');
    UDG = getsolution(app_tdep.buildpath + '/dataout' + "/out_t1",dmd,master.npe);
    mesh.porder = app.porder;
    figure(1); clf; scaplot(mesh, eulereval(UDG, 'M',1.4,0.2)); xlim([-0.1 1.2]); ylim([-0.5 0.5]);
    % [UDG,UH,~,it,alfa,duh] = hdg_solve_dirk(master,mesh,app_dirk,UDG,UH,[],time,dt,nstage,torder);    
    % [UDG,UH] = hdg_solve_dirk(master,mesh,app_dirk,UDG,UH,[],time,dt(itime),nstage,torder);    
    rho_old = sum(UDG_old(:,1:ns,:),2);
    rhoe_old = (UDG_old(:,ns+3,:) - 0.5 * (UDG_old(:,ns+1,:).^2 + UDG_old(:,ns+2,:).^2) ./ rho_old);

    rho_new = sum(UDG(:,1:ns,:),2);
    rhoe_new = (UDG(:,ns+3,:) - 0.5 * (UDG(:,ns+1,:).^2 + UDG(:,ns+2,:).^2) ./ rho_new);

    drho = norm(rho_new(:) - rho_old(:)) / norm(rho_old(:));
    drhoe = norm(rhoe_new(:) - rhoe_old(:)) / norm(rhoe_old(:));
    disp("(dr, dre) = " +string(drho) +", " +string(drhoe) +")");
    time = time + dt;

    alfa = 1;
    % alfa = readbin(app_tdep.buildpath + '/dataout/outalpha.bin');
    if alfa == 1 && drho < 0.1 && drhoe < 0.1
        dt = dt*2;
        disp("Increasing time step: dt = " + string(dt))
    elseif alfa < 0.1 || drho > 1.0 ||drhoe >1.0
        time = time - dt;
        UDG= UDG_old;
        UH = UH_old;
        dt = dt/10;
        disp("Decreasing time step: dt = " + string(dt))
        revert_flag = 1;
    end
    if sum(drho+drhoe) < 1e-5 && ~revert_flag
        % TODO:SAVEALPHA and read it
        disp("Unsteady converged...running steady:")
        mesh.udg = UDG;
        [app,mesh,master,dmd] = preprocessing(app,mesh);
    
        runcode(app, 1); % run C++ code

        % duh_steady = hdg_steady_resid(master,mesh,app_tdep,UDG,UH,[]);
        % disp("Steady residual: " + string(duh_steady));
        % if duh_steady < 1e-7
        disp("Converged")
        break
        % end
    end
    app_tdep.dt = dt*ones(1,1);   % time step sizes
    % app_tdep.read_uh = 1;

    % if it == 1
    %     disp("PTC Converged")
    %     break
    % end 

    % time = time + dt(itime); 
    % save(savedir + "timestep_tmp_new.mat", "UDG", "UH","a","dt","itime");
    % plot_stagline_visc
end
% save(savedir + "timestep_conv_new.mat", "UDG", "UH","a","dt","itime");
end