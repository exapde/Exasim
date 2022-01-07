                     % solution components: u = (u1,...,u_ncu)
xdg = cell(3,1);     % spatial coordinates: xdg = (x_1,..., x_d)
udg = cell(3,1);     % solution and flux:   udg = (u, q) with  q = - nabla u 
odg = cell(3,1);     % other fields:        odg = (odg_1,...,odg_nco)  
uhg = cell(3,1);     % solution trace:      uhg = (uh_1,...,uh_ncu) 
nlg = cell(3,1);     % normal vectors:      nlg = (nl_1,...,nl_d)  
tau = cell(3,1);     % stabilization:       tau = (tau_1,...,tau_ntau)
uinf = cell(3,1);    % boundary data:       uinf = (uinf_1,...,uinf_nuinf)
param = cell(3,1);   % physical parameters: param = (param_1,...,param_npram)

uhg0 = cell(3,1);   
uhg1 = cell(3,1);   
uhg2 = cell(3,1);   
uhg3 = cell(3,1);   
udg0 = cell(3,1);   
udg1 = cell(3,1);   
udg2 = cell(3,1);   
udg3 = cell(3,1);   
odg1 = cell(3,1);
odg2 = cell(3,1);

syms time;  % time 
syms dt;    % timestep size

flux = cell(3,1);    % flux matrix:    flux = [flux_1,..,flux_d], 
                     %                 flux_i = (flux_i1,...,flux_incu), i = 1,...,d 
source = cell(3,1);  % source term:    source = (source_1,..,source_ncu), 
tdfunc = cell(3,1);  % time derivative function:  tdfunc = (tdfunc_1,..,tdfunc_ncu), 
ubou = cell(3,1);    % boundary trace: ubou = (ubou_1,...,ubou_nubou)
fbou = cell(3,1);    % boundary flux:  fbou = (fbou_1,...,fbou_nfbou)
uhat = cell(3,1);    % interior trace: uhat = (uhat_1,...,uhat_ncu)
fhat = cell(3,1);    % interior flux:  fhat = (fhat_1,...,fhat_ncu)
avfd = cell(3,1);    % artificial viscosity
stab = cell(3,1);    % stabilization function
uboutdep = cell(3,1); % time-dependent boundary trace

nd = app.nd;
for d = nd:nd % for each spatial dimension
    ncx = app.ncx;
    xdg{d} = sym('xdg',[ncx 1]); % xdg = (x_1,...,x_d)
    
    ncu = app.ncu;
    ncq = app.ncq; 
    udg{d} = sym('udg',[ncu+ncq 1]); % udg = (u, q)
    udg0{d} = sym('udg0',[ncu+ncq 1]); % udg = (u, q)
    udg1{d} = sym('udg1',[ncu+ncq 1]); % udg = (u, q)
    udg2{d} = sym('udg2',[ncu+ncq 1]); % udg = (u, q)
    udg3{d} = sym('udg3',[ncu+ncq 1]); % udg = (u, q)
    
    nco = app.nco; % no other fields
    odg{d} = sym('odg',[nco 1]);
    odg1{d} = sym('odg1',[nco 1]);
    odg2{d} = sym('odg2',[nco 1]);
    
    ncu = app.ncu;
    uhg{d} = sym('uhg',[ncu 1]); % uhg = (uh_1,...,uh_ncu)
    uhg0{d} = sym('uhg0',[ncu 1]); % uhg = (uh_1,...,uh_ncu)
    uhg1{d} = sym('uhg1',[ncu 1]); % uhg = (uh_1,...,uh_ncu)
    uhg2{d} = sym('uhg2',[ncu 1]); % uhg = (uh_1,...,uh_ncu)
    uhg3{d} = sym('uhg3',[ncu 1]); % uhg = (uh_1,...,uh_ncu)
    
    nlg{d} = sym('nlg',[d 1]); % nlg = (nl_1,...,nl_d)
    
    ntau = length(app.tau);
    tau{d} = sym('tau',[ntau 1]); % tau = (tau_1,...,tau_ntau)
    
    nuinf = length(app.uinf);
    uinf{d} = sym('uinf',[nuinf 1]); % uinf = (uinf_1,...,uinf_nuinf)
        
    nparam = length(app.physicsparam);
    param{d} = sym('param',[nparam 1]); % param = (param_1,...,param_nparam)        
end


