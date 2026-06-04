function avField = getavfield2dchem(y,u,q,hm,~,avcoeff,porder)
%GETAVFIELD2DCHEM Pointwise AV sensor for axisymmetric chemistry model
%
% State:
%   u = [rho_1,...,rho_ns, rho*u_z, rho*u_r, rhoE]
%
% Exasim convention:
%   q = -grad(u)

% number of species / channels
ns  = numel(u) - 3;
nch = ns + 3;

% regularization parameters
alpha = 1.0e3;
rmin  = 1.0e-3;
ymin  = 1.0e-8;

% conservative variables
rho_i = u(1:ns);
rho   = sum(rho_i);
rhou  = u(ns+1);
rhov  = u(ns+2);

% Exasim convention: q = -grad(u)
drho_dz_i = -q(1:ns);
drhou_dz  = -q(ns+1);
drhov_dz  = -q(ns+2);

drho_dr_i = -q(nch+1:nch+ns);
drhou_dr  = -q(nch+ns+1);
drhov_dr  = -q(nch+ns+2);

drho_dz = sum(drho_dz_i);
drho_dr = sum(drho_dr_i);

% regularize rho
rho = rmin + lmax(rho-rmin,alpha);
rhoinv = 1.0 / rho;

% primitive velocities
uz = rhou * rhoinv;
ur = rhov * rhoinv;

% velocity gradients
duz_dz = (drhou_dz - drho_dz * uz) * rhoinv;
dur_dz = (drhov_dz - drho_dz * ur) * rhoinv;

duz_dr = (drhou_dr - drho_dr * uz) * rhoinv;
dur_dr = (drhov_dr - drho_dr * ur) * rhoinv;

% regularize radial coordinate near axis
yreg = ymin + lmax(y-ymin,alpha);

% axisymmetric divergence
divu = duz_dz + dur_dr + ur / yreg;

% compression sensor
comp = -divu;

% limit compression
sigm = 100.0;
comp = limiting(comp,0,sigm,alpha,0);

avField = drhov_dr;

% Optional Ducros-type shear suppression:
% vort_theta = duz_dr - dur_dz;
% vort = sqrt(vort_theta * vort_theta);
% vort = limiting(vort,0,sigm,alpha,0);
% DucrosRatio = comp*comp / (comp*comp + vort*vort + 1.0e-16);

% DucrosRatio = 1.0;
% c_star = 1.0;
% 
% sb = sqrt(hm./porder) * (comp/c_star) * DucrosRatio;
% avField = avcoeff * limiting(sb,0,4,alpha,0.1);
% 
% sb = sqrt(hm./porder) * (comp/c_star) * DucrosRatio;
% 
% % Baseline AV — unchanged everywhere
% av_base = avcoeff * limiting(sb, 0, 4, alpha, 0.1);
% 
% % Extra AV near the nose, fades away quickly
% y_nose   = 0.08;
% w_nose   = 1.0 / (1.0 + (y / y_nose)^2);
% av_nose  = avcoeff * limiting(sb, 0, 4, alpha, 0.1) * w_nose*7.5;
% 
% % Sum: full AV everywhere + bonus near axis
% avField = av_base + av_nose;


end


% function avField = getavfield2dchem(y,u,q,hm,~,avcoeff,porder)
% %GETAVFIELD2DCHEM Pointwise AV sensor for axisymmetric chemistry model
% %
% % State:
% %   u = [rho_1,...,rho_ns, rho*u_z, rho*u_r, rhoE]
% %
% % Exasim convention:
% %   q = -grad(u)
% 
% % number of species / channels
% ns  = numel(u) - 3;
% nch = ns + 3;
% 
% % regularization parameters
% alpha = 1.0e3;
% rmin  = 1.0e-3;
% ymin  = 1.0e-8;
% 
% % conservative variables
% rho_i = u(1:ns);
% rhou  = u(ns+1);
% rhov  = u(ns+2);
% 
% % Reflect negative species densities through zero — abs() is
% % purely arithmetic, no branching, safe for symbolic codegen.
% rho = sum(abs(rho_i));
% 
% % Exasim convention: q = -grad(u)
% drho_dz_i = -q(1:ns);
% drhou_dz  = -q(ns+1);
% drhov_dz  = -q(ns+2);
% drho_dr_i = -q(nch+1:nch+ns);
% drhou_dr  = -q(nch+ns+1);
% drhov_dr  = -q(nch+ns+2);
% 
% % Gradients left as-is — no logical indexing, no masking.
% drho_dz = sum(drho_dz_i);
% drho_dr = sum(drho_dr_i);
% 
% % regularize rho
% rho = rmin + lmax(rho-rmin,alpha);
% rhoinv = 1.0 / rho;
% 
% % primitive velocities
% uz = rhou * rhoinv;
% ur = rhov * rhoinv;
% 
% % velocity gradients
% duz_dz = (drhou_dz - drho_dz * uz) * rhoinv;
% dur_dz = (drhov_dz - drho_dz * ur) * rhoinv;
% duz_dr = (drhou_dr - drho_dr * uz) * rhoinv;
% dur_dr = (drhov_dr - drho_dr * ur) * rhoinv;
% 
% % regularize radial coordinate near axis
% yreg = ymin + lmax(y-ymin,alpha);
% 
% % axisymmetric divergence
% divu = duz_dz + dur_dr + ur / yreg;
% 
% % compression sensor
% comp = -divu;
% 
% % limit compression
% sigm = 100.0;
% comp = limiting(comp,0,sigm,alpha,0);
% 
% DucrosRatio = 1.0;
% c_star = 1.0;
% 
% sb = sqrt(hm./porder) * (comp/c_star) * DucrosRatio;
% av_base = avcoeff * limiting(sb, 0, 6, alpha, 0.1);
% 
% avField = av_base;
% end
