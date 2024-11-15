% syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21 u22 u23 u24
% syms w1 w2
% syms dw1_du1 dw1_du2 dw1_du3 dw1_du4 dw1_du5 dw1_du6 dw1_du7 dw1_du8
clear
syms x1 x2 x3 av
syms time
syms param1 param2 param3 param4 param5 param6 param7 param8 param9 param10 param11 param12 param13
syms zero one 

param = [param1 param2 param3 param4 param5 param6 param7 param8 param9 param10 param11 param12 param13];
nd=2;

% udg = [u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21 u22 u23 u24]; 
udg = sym("u",[1 24]); % (ns+nd+1)*(nd+1)
wdg = sym("w",[1 2]); % needs to be at least 2
dw_du = sym("dw_du",[2 8]);  
% wdg = [w1 w2];
% dw_du = [dw1_du1 dw1_du2 dw1_du3 dw1_du4 dw1_du5 dw1_du6 dw1_du7 dw1_du8];
pg = [x1 x2 x3];       

ns = 5;
ncu = ns + nd + 1;
nc = length(udg);
nch = ncu;

% Nondimensional params
rho_scale   = param(1);
u_scale     = param(2);
rhoe_scale  = param(3);
T_scale     = param(4);
mu_scale    = param(5);
kappa_scale = param(6);
cp_scale    = param(7);
L_scale     = param(8);
% Ec          = param(9);
Ec = 1.0;
Pr          = param(10);
Re          = param(11);
beta = 0;

rho_i = udg(1:ns);
rho   = sum(rho_i);
rhou  = udg(ns+1);
rhov  = udg(ns+2);
rhoE  = udg(ns+3);

drho_dx_i = -udg(nch + (1:ns));
drhou_dx  = -udg(nch + ns+1);
drhov_dx  = -udg(nch + ns+2);
drhoE_dx  = -udg(nch + ns+2+1);
drho_dy_i = -udg(nch + (nch+1:nch+ns));
drhou_dy  = -udg(nch + nch+ns+1);
drhov_dy  = -udg(nch + nch+ns+2);
drhoE_dy  = -udg(nch + nch+ns+2+1);
av = pg(3);

% alphaClip = 1e12;
% rmin = 0;
% Conservative Variables
% for i = 1:ns
%     rho_i(i) = rmin + lmax(udg(i)-rmin,alphaClip); %subspecies density
%     rho = rho + rho_i(i); %total mixture density
% end

rho_inv = 1.0 ./ rho;
% uv = rhou .* rho_inv; %velocity
% vv = rhov .* rho_inv;
%energy
% uTu2   = 0.5*(uv.*uv+vv.*vv);
% Y_i = rho_i / 

% i_P = 1;
% i_T = 2;
% i_dT_drhoi = 3:ns+2;
% i_dT_drhoe = ns+3;

% p = wdg(i_P);
% T = wdg(i_T); 
% dT_drho_i = wdg(i_dT_drhoi);
% dT_drhoe = wdg(i_dT_drhoe);
i_T = 1;
i_P = 2;
T = wdg(i_T);
p = wdg(i_P);
T_dim = T*T_scale;
p_dim = p*rhoe_scale;

C1 = 1.458e-6;
C2 = 110.4;

% mu_S = C1 * T.^(3/2) ./ (C2 + T);
[species_thermo_structs, Mw, ~] = thermodynamicsModels();
[blottner_structs, gupta_structs, gupta_mu_structs, gupta_kappa_structs] = transport();

% Mw = Mw';

% mu_d_dim = C1 * T_dim.^(3/2) ./ (C2 + T_dim);
% kappa_dim = mu_d_dim * cp_scale / Pr;

rho_i_dim = rho_i*rho_scale;
Y = Y_i(rho_i_dim);
denom = sum(rho_i_dim) * mixtureFrozenCvMass(T_dim, Mw, Y, species_thermo_structs);
e_i = getEnergiesMass(T_dim, Mw, species_thermo_structs);
dT_drho_i_dim = -e_i(:)' ./ denom;
dT_drhoe_dim = 1.0 / denom;

dT_drho_i = dT_drho_i_dim / T_scale * rho_scale;
dT_drhoe = dT_drhoe_dim / T_scale * rhoe_scale;

% Some useful derived quantities 
rho = sum(rho_i);
drho_dx = sum(drho_dx_i);
drho_dy = sum(drho_dy_i);
rho_inv = 1.0 / rho;
X = X_i(rho_i_dim,Mw);
%     drho_dx_inv = 1.0 / drho_dx;
%     drho_dy_inv = 1.0 / drho_dy
uv = rhou * rho_inv; %velocity
vv = rhov * rho_inv;
E = rhoE .* rho_inv; 
du_dx = (drhou_dx - drho_dx*uv)*rho_inv;
dv_dx = (drhov_dx - drho_dx*vv)*rho_inv;
du_dy = (drhou_dy - drho_dy*uv)*rho_inv;
dv_dy = (drhov_dy - drho_dy*vv)*rho_inv;
% re = rhoE - rho * (0.5 * uTu); 
uTu2      = 0.5 * (uv * uv + vv * vv);
duTu2_dx  = uv * du_dx + vv * dv_dx; 
duTu2_dy  = uv * du_dy + vv * dv_dy;
%%% TODO: ECKERT NUMBER HERE? 
dre_drho  = -uTu2;
dre_duTu2 = -rho;
dre_drhoE = 1.0;
dre_dx    = dre_drho * drho_dx + dre_duTu2 * duTu2_dx + dre_drhoE * drhoE_dx;
dre_dy    = dre_drho * drho_dy + dre_duTu2 * duTu2_dy + dre_drhoE * drhoE_dy;
dT_dx     = sum(dT_drho_i .* drho_dx_i) +  dT_drhoe * dre_dx;
dT_dy     = sum(dT_drho_i .* drho_dy_i) +  dT_drhoe * dre_dy;

H = E + p.*rho_inv; %enthalpy

D_vec = averageDiffusionCoeffs(T_dim, X, Y, Mw, p_dim, gupta_structs);
h_vec = getEnthalpiesMass(T_dim, Mw, species_thermo_structs);

mu_i = speciesViscosities(T_dim, gupta_mu_structs);
% lambda_i = speciesConductivities(T, gupta_kappa_structs);
phi_i = euckenPhi(mu_i, Mw, X);
mu_d_dim = wilkeMixture(mu_i, X, phi_i);
lambda_i = mu_d_dim * 3/2 .* getCpsMass(T_dim, Mw, species_thermo_structs);
kappa_dim = wilkeMixture(lambda_i, X, phi_i);

h_scale = u_scale^2;
D_scale = u_scale;

mu_d = mu_d_dim / mu_scale;
kappa = kappa_dim / kappa_scale;
D_vec = D_vec' / D_scale;
h_vec = h_vec' / h_scale;

%%%%%%%% Calculation of J_i
%     Y_i = rho_i ./ rho;
dY_dx_i = (drho_dx_i * rho - rho_i * drho_dx) * rho_inv * rho_inv;
dY_dy_i = (drho_dy_i * rho - rho_i * drho_dy) * rho_inv * rho_inv;

J_i_x = -rho * D_vec .* dY_dx_i + rho_i .* sum(D_vec .* dY_dx_i);
J_i_y = -rho * D_vec .* dY_dy_i + rho_i .* sum(D_vec .* dY_dy_i);

%%%%%%%% Stress tensor tau
% txx = viscshear * 4.0/3.0 * du_dx + viscbulk * du_dx; %TODO: check sign
txx = mu_d * 2.0/3.0 * (2 * du_dx - dv_dy) / Re + beta * (du_dx + dv_dy);
txy = mu_d * (du_dy + dv_dx) / Re;
tyy = mu_d * 2.0/3.0 * (2 * dv_dy - du_dx) / Re + beta * (du_dx + dv_dy);


% Fluxes
fi = sym(zeros(nch,2));
fv = sym(zeros(nch,2));

for i = 1:ns
    fi(i,1) = rho_i(:,i) .* uv     - av.*drho_dx_i(:,i);
end
fi(ns + 1,1) = rhou .* uv + p    - av.*drhou_dx;
fi(ns + 2,1) = rhov .* uv        - av.*drhov_dx;
fi(ns + 3,1) = rhou .* H  - av.*drhoE_dx;

for i = 1:ns
    fi(i,2) = rho_i(:,i) .* vv     - av.*drho_dy_i(:,i);
end
fi(ns + 1,2) = rhou .* vv        - av.*drhou_dy;
fi(ns + 2,2) = rhov .* vv + p    - av.*drhov_dy;
fi(ns + 3,2) = rhov .* H         - av.*drhoE_dy;
    
%%% TODO: MANUALLY TRIPLE CHECK AV FIELD SIGN 
f2 = fi(:);
for n=2
    if n==1
        f = f1;
        filename1 = ['fluxvisc' num2str(nd) 'd' '.m'];
    elseif n==2
        f = f2;
        filename1 = ['flux' num2str(nd) 'd2_mpp_manual' '.m'];
    end
    
    %%% compute Jacobian
    jac_f = jacobian(f,udg);
    jac_w = jacobian(f, wdg);
    %%% And patch with vector zero or one to have the right sizes
    for ii = 1:size(f,1)
        for jj = 1:size(f,2)
            temp = ismember(symvar(f(ii,jj)),udg);
            if f(ii,jj)==0, f(ii,jj) = zero; %%% No dependency on anything
            elseif isempty(temp) || sum(temp)==0, f(ii,jj) = (f(ii,jj))*one; %%% No dependency on state vars
            end
        end
    end

    for ii = 1:size(jac_f,1)
        for jj = 1:size(jac_f,2)
            temp = ismember(symvar(jac_f(ii,jj)),udg);
            if jac_f(ii,jj)==0, jac_f(ii,jj) = zero; %%% No dependency on anything
            elseif isempty(temp) || sum(temp)==0, jac_f(ii,jj) = (jac_f(ii,jj))*one; %%% No dependency on state vars
            end
        end 
    end    
    
    %matlabFunction(f(1),'file','tmp.m','vars',{pg,udg,param,time,[zero one]},'outputs', {'f'});
    
    % generate a temporary matlab file
    jac_f(:,1:ncu) = jac_f(:,1:ncu) + jac_w * dw_du;
    matlabFunction(f(:),jac_f(:),'file','tmp.m','vars',{pg,udg,wdg,dw_du,param,time,[zero one]},'outputs', {'f','f_udg'});

    % open the file and modify it
    fid = fopen('tmp.m','r');
    gid = fopen(filename1,'wt');

    tline = fgetl(fid); 
    i=1;       
    while ischar(tline)        
        str = strrep(tline, 'tmp', strrep(filename1,'.m',''));
        str = strrep(str, 'TMP', upper(strrep(filename1,'.m','')));
        str = strrep(str, 'in1', 'pg');
        str = strrep(str, 'IN1', 'PG');        
        str = strrep(str, 'in2', 'udg');
        str = strrep(str, 'IN2', 'UDG'); 
        str = strrep(str, 'in3', 'wdg');
        str = strrep(str, 'IN3', 'WDG');     
        str = strrep(str, 'in4', 'dw_du');
        str = strrep(str, 'IN4', 'DW_DU');
        str = strrep(str, 'in5', 'param');
        str = strrep(str, 'IN5', 'PARAM');     
        str = strrep(str, ',in7)', ')');                
        str = strrep(str, ',IN7)', ')');    
        str = strrep(str, 'param(:,1)', 'param{1}'); 
        str = strrep(str, 'param(:,2)', 'param{2}'); 
        str = strrep(str, 'param(:,3)', 'param{3}'); 
        str = strrep(str, 'param(:,4)', 'param{4}'); 
        str = strrep(str, 'param(:,5)', 'param{5}'); 
        str = strrep(str, 'param(:,6)', 'param{6}'); 
        str = strrep(str, 'param(:,7)', 'param{7}'); 
        str = strrep(str, 'param(:,8)', 'param{8}'); 
        str = strrep(str, 'param(:,9)', 'param{9}');     
        str = strrep(str, 'param(:,10)', 'param{10}');  
        str = strrep(str, 'param(:,11)', 'param{11}');     
        str = strrep(str, 'in7(:,1)', 'zeros(ng,1)');
        str = strrep(str, 'in7(:,2)', 'ones(ng,1)');        
        if i==7
            str = '[ng,nc] = size(udg);';
            fprintf(gid, '%s\n', str);                  
            str = ['nch = ' num2str(nch) ';'];
            fprintf(gid, '%s\n', str);                  
            str = ['nd = ' num2str(nd) ';'];
        end
        fprintf(gid, '%s\n', str);                  
        tline = fgetl(fid);        
        i=i+1;
        %disp(str)
    end

    str = 'f = reshape(f,ng,nch,nd);';
    fprintf(gid, '%s\n', str);                  
    str = 'f_udg = reshape(f_udg,ng,nch,nd,nc);';
    fprintf(gid, '%s\n', str);                  

    fclose(fid);
    fclose(gid);
    delete('tmp.m');
end

