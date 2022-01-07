syminit; % initialize symbolic variables
f = getflux4_3d(udg{d},odg{d},param{d});
% define flux as a function of xdg, udg, odg, param, time
flux{d} = reshape(f,ncu,d);
% define source as a function of xdg, udg, odg, param, time
source{d} = [0.0; 0.0; 0.0; 0.0; sym(0.0)];
% define tdfunc as a function of xdg, udg, odg, param, time
tdfunc{d} = [1.0; 1.0; 1.0; 1.0; sym(1.0)];
u = udg{d}(1:ncu);
ui = uinf{d}(1:ncu);
An = getan_3d(nlg{d},u,param{d}(1),0);
ub = 0.5*((u+ui) + An*(u-ui));
uw = u;
uw(2:4) = 0;
% define ubou as a function of xdg, udg, odg, uhg, nlg, tau, uinf, param, time
ubou{d} = [uw ub];
% define fbou as a function of xdg, udg, odg, uhg, nlg, tau, uinf, param, time
f = getflux_3d([uhg{d}; udg{d}(ncu+1:end)],param{d});
f = reshape(f,ncu,d);
fx = f(:,1);
fy = f(:,2);
fz = f(:,3);
nx = nlg{d}(1);
ny = nlg{d}(2);
nz = nlg{d}(3);
fb = fx*nx + fy*ny + fz*nz + tau{d}*(u-uhg{d});
fw = fb;
fw(end) = 0.0;
fbou{d} = [fw fb];
% generate AV field
avfd{d} = getavfield2_3d(udg{d},odg{d},param{d},app.porder);
% produce executable programs
gencode;
