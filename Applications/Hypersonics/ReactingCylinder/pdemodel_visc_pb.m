 
function pde = pdemodel
    pde.mass = @mass;
    pde.flux = @flux;
    pde.source = @source;
    pde.fbou = @fbou;
    pde.ubou = @ubou;
    pde.initu = @initu;
    pde.output = @output;
    pde.avfield = @avfield;
end

function f = avfield(u, q, w, v, x, t, mu, eta)
    nspecies = 5;
    nd = 2;
    f = avfieldPhysical2d_new(u, q, w, v, x, t, mu, eta, nspecies, nd);
end

function m = mass(u, q, w, v, x, t, mu, eta)
    nspecies = 5;
    ndim = 2;
    m = sym(ones(nspecies + ndim + 1, 1));
end

function f = flux(u, q, w, v, x, t, mu, eta)
    nspecies = 5;
    nd = 2;

    fi = fluxinviscid(u, q, w, v, x, t, mu, eta, nspecies, nd);
    fv = fluxviscousPb(u, q, w, v, x, t, mu, eta, nspecies, nd);
    f = fi - fv;
end

function fi = fluxinviscid(u, q, w, v, x, t, mu, eta, nspecies, ndim)
    nenergy = 1;
    % ndim = 2;
    fi = sym(zeros(nspecies+ndim+1,2));

    rho_i = sym(zeros(nspecies,1));
    rho = sym(0);

    rmin = 0.0;
    alpha = mu(22);

    % Conservative Variables
    for ispecies = 1:nspecies
        rho_i(ispecies) = rmin + lmax(u(ispecies)-rmin,alpha); %subspecies density
        rho = rho + rho_i(ispecies); %total mixture density
    end

    rhou = u(nspecies+1);
    if ndim == 2
        rhov = u(nspecies+2);
    end
    rhoE = u(nspecies+ndim+1);
    if nenergy == 2
        rhoev = u(nspecies+ndim+nenergy);
    end

    rhoinv = 1.0 / rho;
    uv = rhou * rhoinv; %velocity
    vv = rhov * rhoinv;
    E = rhoE * rhoinv; %energy

    p = eta(1); % MANUALLY MODIFIED
    H = E + p*rhoinv; %enthalpy

    % Fluxes
    for ispecies = 1:nspecies
        fi(ispecies,1) = rho_i(ispecies) * uv;
    end
    fi(nspecies + 1,1) = rhou * uv + p;
    fi(nspecies + 2,1) = rhov * uv;
    fi(nspecies + ndim + 1,1) = rhou * H;

    for ispecies = 1:nspecies
        fi(ispecies,2) = rho_i(ispecies) * vv;
    end
    fi(nspecies + 1,2) = rhou * vv;
    fi(nspecies + 2,2) = rhov * vv + p;
    fi(nspecies + ndim + 1,2) = rhov * H;

    % fi = [ru, ru*uv+p, rv*uv, ru*h, ...
            % rv, ru*vv, rv*vv+p, rv*h];
    % fi = reshape(fi,[4,2]);    
end



function fv = fluxviscousPb(u, q, w, v, x, t, mu, eta, nspecies, ndim)
    fv = sym(zeros(nspecies+ndim+1,2));
    rmin = 0.0;
    alpha = mu(22);
    nenergy=1;
    ncu = nspecies + ndim + nenergy;

    Pr = mu(19);
    Re = mu(20);

    avb = v(1);
    avk = mu(21) * avb;
    avs = 0;
    avk_b_coeff = eta(3*nspecies+5);

    % Data layout
    %   [D_1, ..., D_ns, h_1, ..., h_ns, mu, kappa, dTdri, dTdrhou, dTdrhoE]
    p = eta(1);
    D_vec = eta(2:nspecies+1); 
    h_vec = eta(nspecies+2:2*nspecies+1);
    mu_d = eta(2*nspecies+2);
    kappa = eta(2*nspecies+3);
    dT_drho_i = eta(2*nspecies+4:3*nspecies+3);
    dT_drhoe = eta(3*nspecies+4);
    beta = 0.0;

    rho_i = rmin + lmax(u(1:nspecies) - rmin, alpha);
    dr = dlmax(u(1:nspecies)-rmin,alpha); % Regularize derivative 

    rhou = u(nspecies+1);
    rhov = u(nspecies+2);
    rhoE = u(nspecies+ndim+1);
    
    drho_dx_i = -q(1:nspecies) .* dr;
    drhou_dx  = -q(nspecies+1);
    drhov_dx  = -q(nspecies+2);
    drhoE_dx  = -q(nspecies+ndim+1);
    drho_dy_i = -q(ncu+1:ncu+nspecies) .* dr;
    drhou_dy  = -q(ncu+nspecies+1);
    drhov_dy  = -q(ncu+nspecies+2);
    drhoE_dy  = -q(ncu+nspecies+ndim+1);

    % Some useful derived quantities 
    rho = sum(rho_i);
    drho_dx = sum(drho_dx_i);
    drho_dy = sum(drho_dy_i);
    rho_inv = 1.0 / rho;

    uv = rhou * rho_inv; %velocity
    vv = rhov * rho_inv;
    du_dx = (drhou_dx - drho_dx*uv)*rho_inv;
    dv_dx = (drhov_dx - drho_dx*vv)*rho_inv;
    du_dy = (drhou_dy - drho_dy*uv)*rho_inv;
    dv_dy = (drhov_dy - drho_dy*vv)*rho_inv;

    % re = rhoE - rho * (0.5 * uTu); 
    uTu2      = 0.5 * (uv * uv + vv * vv);
    duTu2_dx  = uv * du_dx + vv * dv_dx; 
    duTu2_dy  = uv * du_dy + vv * dv_dy;
    dre_drho  = -uTu2;
    dre_duTu2 = -rho;
    dre_drhoE = 1.0;
    dre_dx    = dre_drho * drho_dx + dre_duTu2 * duTu2_dx + dre_drhoE * drhoE_dx;
    dre_dy    = dre_drho * drho_dy + dre_duTu2 * duTu2_dy + dre_drhoE * drhoE_dy;
    dT_dx     = sum(dT_drho_i .* drho_dx_i) +  dT_drhoe * dre_dx;
    dT_dy     = sum(dT_drho_i .* drho_dy_i) +  dT_drhoe * dre_dy;

    %%%%%%%% Modify appropriate coefficients 
    kappa = kappa / (Re * Pr) + avk * avk_b_coeff; % thermal conductivity
    mu_d = mu_d / Re + avs;   % dynamic viscosity
    beta = 0 + avb;  % bulk viscosity 
    D_vec = D_vec + 0; % Species diffusion velocities

    %%%%%%%% Calculation of J_i
%     Y_i = rho_i ./ rho;
    dY_dx_i = (drho_dx_i * rho - rho_i * drho_dx) * rho_inv * rho_inv;
    dY_dy_i = (drho_dy_i * rho - rho_i * drho_dy) * rho_inv * rho_inv;
 
    J_i_x = -rho * D_vec .* dY_dx_i + rho_i .* sum(D_vec .* dY_dx_i);
    J_i_y = -rho * D_vec .* dY_dy_i + rho_i .* sum(D_vec .* dY_dy_i);

    %%%%%%%% Stress tensor tau
    % txx = viscshear * 4.0/3.0 * du_dx + viscbulk * du_dx; %TODO: check sign
    txx = mu_d * 2.0/3.0 * (2 * du_dx - dv_dy) + beta * (du_dx + dv_dy);
    txy = mu_d * (du_dy + dv_dx);
    tyy = mu_d * 2.0/3.0 * (2 * dv_dy - du_dx) + beta * (du_dx + dv_dy);
    %%%%%%%% Final viscous fluxes
    for ispecies = 1:nspecies
        fv(ispecies,1) = -J_i_x(ispecies);
        fv(ispecies,2) = -J_i_y(ispecies);
    end

    fv(nspecies + 1, 1) = txx;
    fv(nspecies + 2, 1) = txy;
    fv(nspecies + 3,1) = uv * txx + vv * txy - (sum(h_vec.*J_i_x) - kappa * dT_dx);
    fv(nspecies + 1, 2) = txy;
    fv(nspecies + 2, 2) = tyy;
    fv(nspecies + 3,2) = uv * txy + vv * tyy - (sum(h_vec.*J_i_y) - kappa * dT_dy);
end


function s = source(u, q, w, v, x, t, mu, eta)
    
    nspecies = 5;
    rmin = 0.0;
    alpha = mu(22);
    rvec = u(1:nspecies);
    s(1:nspecies) = rmin + lmax(rvec-rmin,alpha);
    s(nspecies+1) = 0.0;
    s(nspecies+2) = 0.0;
    s(nspecies+3) = 0.0;
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    nspecies = 5;
    ndim = 2;
    ub = sym(zeros(nspecies+ndim+1, 2));
%%%%%% C-CODE MANUALLY WRITTEN
%     ub(:,1) = u(:);
%     ub(:,2) = u(:); 

    uinflow = initu(x, mu, eta);

    uoutflow = u;

    uadiabatic = u;
    uadiabatic(nspecies+1:nspecies+ndim) = 0.0;

    ub(:,1) = uinflow;
    ub(:,2) = uoutflow;
    ub(:,3) = uadiabatic;
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    nspecies = 5;
    ndim = 2;
    fb = sym(zeros(nspecies+ndim+1, 3));

    f = flux(u, q, w, v, x, t, mu, eta);
    fn = f(:,1)*n(1) + f(:,2)*n(2) + tau.*(u-uhat);


    fadiabatic = fn;
    fadiabatic(1:nspecies) = 0.0;
    fadiabatic(nspecies+ndim+1) = 0.0;

    fb(:,1) = fn;
    fb(:,2) = fn;
    fb(:,3) = fadiabatic;
end

function u0 = initu(x, mu, eta)
    nspecies = 5;
    rho_scale = eta(1);
    u_scale = eta(2);
    rhoe_scale = eta(3);

%     u0 = mu(1:nspecies+2+1);
    u0(1:nspecies) = mu(1:nspecies) / rho_scale;
    u0(nspecies+1:nspecies+2) = mu(nspecies+1:nspecies+2) / (rho_scale * u_scale);
    u0(nspecies+2+1) = mu(nspecies + 2 + 1) / rhoe_scale;
end

function o = output(u, q, w, v, x, t, mu, eta)
    
    nspecies = 5;

    o = sym(zeros(10,1));
    nspecies = 5;
    alpha = mu(22);
    rmin = 0.0;

    f = sym(zeros(nspecies + 2,1));
    rho_i = sym(zeros(5,1));
    rho = sym(0);

    % Conservative Variables
    for ispecies = 1:nspecies
        rho_i(ispecies) = rmin + lmax(u(ispecies)-rmin,alpha); %subspecies density
        rho = rho + rho_i(ispecies); %total mixture density
    end

    rhou = u(nspecies+1);
    rhoE = u(nspecies+2);

    rhoinv = 1.0 / rho;
    uv = rhou * rhoinv; %velocity
    E = rhoE*rhoinv; %energy

    o(1:8) = u(1:8);
    o(1) = v(1);
%     o = sensor_outputs(u, q, w, v, x, t, mu, eta);
end

% function dout = dlmax(x, alpha)
%     dout = atan(alpha*(x))/pi + (alpha*(x))./(pi*(alpha^2*(x).^2 + 1)) + 1/2;
% end