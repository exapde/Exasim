syminit; % initialize symbolic variables
f = getfluxav3d(udg{d},odg{d},param{d});
% define flux as a function of xdg, udg, odg, param, time
flux{d} = reshape(f,ncu,d);
% define source as a function of xdg, udg, odg, param, time
source{d} = [0.0; 0.0; 0.0; 0.0; sym(0.0)];
% define tdfunc as a function of xdg, udg, odg, param, time
tdfunc{d} = [1.0; 1.0; 1.0; 1.0; sym(1.0)];
u  = udg{d}(1:ncu);
uh = uhg{d}(1:ncu);
ui = uinf{d}(1:ncu);
nx = nlg{d}(1);
ny = nlg{d}(2);
nz = nlg{d}(3);
An = getan3d(nlg{d},u,param{d}(1),0);
% Freestream BC
ub = 0.5*((u+ui) + An*(u-ui));
% Wall BC
uw = u;
uw(2:4) = 0;
% Slip wall BC
usw = u;
usw(2) = u(2) - nx * (u(2)*nx + u(3)*ny + u(4)*nz);
usw(3) = u(3) - ny * (u(2)*nx + u(3)*ny + u(4)*nz);
usw(4) = u(4) - nz * (u(2)*nx + u(3)*ny + u(4)*nz);
% Isothermal Wall
%TisoW = Twall/Tref * Tinf;
p_inf = 1/(param{d}(1)*(param{d}(4)^2));
T_inf = p_inf/(param{d}(1)-1);
TisoW = param{d}(10)/param{d}(9) * T_inf;
utw = u;
utw(2:4) = 0;
utw(5) = u(1)*TisoW;
uo = get_outletBC(u,param{d},nlg{d});
%ut(4) = TisoW - gam*(gam-1.)*Minf*Minf*u(4)/u(1); Pablo's isoth BC
% define ubou as a function of xdg, udg, odg, uhg, nlg, tau, uinf, param, time
%ubou{d} = [ui u ub usw uw]; % [Inflow, Outflow, FreeStream, SlipWall, Wall]
ubou{d} = [ub ui u usw uw utw uo];  % [FreeStream, Inflow, Outflow, SlipWall, Wall, TherWall]
% define fbou as a function of xdg, udg, odg, uhg, nlg, tau, uinf, param, time
f = getfluxav3d([uhg{d}; udg{d}(ncu+1:end)],odg{d},param{d});
%f = getflux([uhg{d}; udg{d}(ncu+1:end)],param{d});
f = reshape(f,ncu,d);
fx = f(:,1);
fy = f(:,2);
fz = f(:,3);
fb = fx*nx + fy*ny + fz*nz +  tau{d}* (u-uhg{d});
% Flux on adiabatic Wall
fw = fb;
fw(1) = 0.0;
fw(end) = 0.0;
% Flux Thermal Wall
%Am = getan3d(nlg{d},u,param{d}(1),1);
%ftw = fx*nx + fy*ny + fz*nz +  Am* (u-uhg{d});
ftw = fb;
ftw(1) = 0.0;
fbou{d} = [fb fb fb fw fw ftw fb];
% generate AV field
avfd{d} = getavfield3d(udg{d},odg{d},param{d},app.porder);
% produce executable programs
gencode;


