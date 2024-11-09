function [w, dw_du] = state2mutation_my_temp(udg, param)

%Input: U fluid state
%output: pressure_mat, dpressure_mat_dU

T_0 = 500;

% Nondimensional params
rho_scale   = param{1};
u_scale     = param{2};
rhoe_scale  = param{3};
T_scale     = param{4};
mu_scale    = param{5};
kappa_scale = param{6};
cp_scale    = param{7};
L_scale     = param{8};
Ec = param{9};
[species_thermo_structs, Mw, RU] = thermodynamicsModels();

ng = size(udg, 1);
ncu = 8;
ns = 5;

i_T = 1;
i_P = 2;


nw = 2; %CHANGE
w = zeros(ng, nw);
dw_du = zeros(ng, nw, ncu);

rho_i = udg(:,1:ns);
rho   = sum(rho_i,2);
rhou  = udg(:,ns+1);
rhov  = udg(:,ns+2);
rhoE  = udg(:,ns+3);

rhoinv = 1.0 ./ rho;
uv = rhou .* rhoinv; %velocity
vv = rhov .* rhoinv;
E = rhoE .* rhoinv; %energy
uTu2   = 0.5*(uv.*uv+vv.*vv);

% TODO: PACK INPUTS
rho_i_dim = rho_i * rho_scale;
rho_dim = sum(rho_i_dim,2);
rhou_dim = rhou * (rho_scale * u_scale);
rhov_dim = rhov * (rho_scale * u_scale);
rhoE_dim = rhoE * rhoe_scale;
rhoe_dim = (rhoE_dim - 0.5 * (rhou_dim.*rhou_dim + rhov_dim.*rhov_dim) ./ rho_dim);

% TODO: PASS TO MUTATION
[species_thermo_structs, Mw, RU] = thermodynamicsModels();
rho_tilde = rho_i_dim ./ Mw';
alpha = -sum(rho_tilde, 2);
T_2 = zeros(size(udg,1),1);
P_2 = zeros(size(udg,1),1);
dT_dr_2 = zeros(size(udg,1),ns);
dT_dre_2 = zeros(size(udg,1),1);
% T_0 = 1000;
% Uin = [rho_i_dim, rhoe_dim];
% disp(Uin(1,:))
% [P, T,dT_dr,dT_dre] = mppFlux_rre(ng, 2, Uin);

for ig = 1:size(rho_i_dim,1)
    T_0 = 1000;
    iter = 1;
    f0 = f_T(T_0, rho_tilde(ig,:), rhoe_dim(ig,1), alpha(ig), species_thermo_structs);
    while max(abs(f0)) > 1e-9
        cp = cp_n(T_0, species_thermo_structs);
        df_dw = alpha(ig) + dot(rho_tilde(ig,:),cp);
        dT = f0 / df_dw;
        while T_0 - dT < 50
            dT = dT * 0.5; % simple line search
        end
        T_0 = T_0 - dT;
        f0 = f_T(T_0, rho_tilde(ig,:), rhoe_dim(ig,1), alpha(ig), species_thermo_structs);
        iter = iter + 1;
        if iter > 100
            % disp("Over max iters...")
            % disp(f0)
            % disp(T_0)
            break
        end
        if T_0 < 50
            T_0 = 50;
        end
    end
    T_2(ig) = T_0;
    P_2(ig) = pressure_mat(T_0, rho_i_dim(ig,:), Mw');
    Y = Y_i(rho_i_dim(ig,:));
    denom = sum(rho_i_dim(ig,:)) * mixtureFrozenCvMass(T_0, Mw, Y, species_thermo_structs);
    e_i = getEnergiesMass(T_0, Mw, species_thermo_structs);
    dT_dr_2(ig,:) = -e_i(:)' ./ denom;
    dT_dre_2(ig) = 1.0 / denom;
end
% disp(norm(T-T_2))
% figure(1); clf; plot(T); hold on; plot(T_2, '--')
% [~, ~, ~, ~,dT_dr, dT_dre,~, ~,~, ~,~, ~,~, ~, ~, ~,~, ~, P, T] = mppFlux_visc(ng, 1e-6, Uin);\
% [p,T,dT_dr,dT_dre] = mpp
P = P_2; 
T = T_2; 
dT_dr = dT_dr_2;
dT_dre = dT_dre_2;
% CHANGE OF VARIABLES
denom = 1 ./ dT_dre;
dT_drho_i = (uTu2*u_scale^2 ./denom) + dT_dr;
dT_drhou = -uv ./ denom * u_scale;
dT_drhov = -vv ./ denom * u_scale;
dT_drhoE = dT_dre;

Mw = [14.0067...
15.9994...
30.0061...
28.0134...
31.9988]; % mix.speciesMw
RU = 8.314471468617452;

Mw = Mw / 1000.0;
dp_drho_i = (T .* RU./Mw  + pressure_mat(dT_drho_i, rho_i_dim, Mw)) / rhoe_scale * rho_scale;
dp_drhou = pressure_mat(dT_drhou, rho_i_dim, Mw) / rhoe_scale * rho_scale * u_scale;
dp_drhov = pressure_mat(dT_drhov, rho_i_dim, Mw) / rhoe_scale * rho_scale * u_scale;
dp_drhoE = pressure_mat(dT_drhoE, rho_i_dim, Mw) / rhoe_scale * rhoe_scale;

w(:,i_P) = P / rhoe_scale;
w(:,i_T) = T / T_scale;

% pressure_mat derivatives
dw_du(:,i_P,1) = dp_drho_i(:,1);
dw_du(:,i_P,2) = dp_drho_i(:,2);
dw_du(:,i_P,3) = dp_drho_i(:,3);
dw_du(:,i_P,4) = dp_drho_i(:,4);
dw_du(:,i_P,5) = dp_drho_i(:,5);
dw_du(:,i_P,6) = dp_drhou;
dw_du(:,i_P,7) = dp_drhov;
dw_du(:,i_P,8) = dp_drhoE;

% Temperature derivatives
dw_du(:,i_T,1:5) = dT_drho_i / T_scale * rho_scale;
dw_du(:,i_T,6) = dT_drhou / T_scale *  rho_scale * u_scale;
dw_du(:,i_T,7) = dT_drhov / T_scale * rho_scale * u_scale;
dw_du(:,i_T,8) = dT_drhoE / T_scale * rhoe_scale;


end