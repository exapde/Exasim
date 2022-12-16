 
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
    f = avfieldLaplacian(u, q, w, v, x, t, mu, eta, nspecies, nd);
end

function avField = avfieldLaplacian(u, q, w, v, x, t, mu, eta, nspecies, nd)
    % Mutation outputs
    p = eta(1);
    gam = eta(2);
    gam1 = gam - 1;
    
    % artificial viscosity
    porder = mu(15);

    % regularization parameters for the bulk viscosity
    kb = mu(16);
    sb0   = mu(17); %0.02 
    sbmin = 0.0;
    sbmax = mu(18);% / sqrt(gam*gam - 1.0); %2.5 TODO: unsure that this should be changing 

    % regularization parameters
    alpha = 1.0e12;
    rmin = 0.0;
    Hmin = 1.0e-4;

    % mesh size
    hm = v(2);

    % Get base variables
    rvec = u(1:nspecies);
    ru = u(nspecies+1);
    rv = u(nspecies+nd);
    rE = u(nspecies+nd+1);

    rx = sum(-q(1:nspecies));
    rux = -q(nspecies+1);
    rvx = -q(nspecies+nd);

    % Regularization 
    rvec = rmin + lmax(rvec-rmin,alpha); % need to double check regularizatino here
    
    % Simple derived working quantities
    r = sum(rvec);
    r1 = 1./r;
    uv = ru.*r1;
    vv = rv.*r1;
    E = rE * r1; %energy
    H = E + p*r1; %enthalpy
    H = Hmin + lmax(H - Hmin, alpha);

    % Critical speed of Sound
    c_star = sqrt((2.*gam1.*H) ./ (gam+1)); %TODO: this is fine for 1 temp but not 2

    % Computing derivatives for the sensors
    ux = (rux - rx.*uv).*r1;
    vx = (rvx - rx.*vv).*r1;
    div_v = (ux + vx); %TODO: check signs. Actually this probably is important, we want to make sure it's applied in negative v. pos dilitation
                  %      pretty sure my signs are okay. There is something strange about the code I was given. Truly don't understand the signs 
    % limit  divergence and vorticity
    sigm = 1e4;
    div_v = limiting(div_v,-sigm,sigm,alpha,-sigm);

    % % Dilatation Sensor sb
    sb = - (hm./porder) .* (div_v./c_star);
    sb = limiting(sb,sbmin,sbmax,alpha,sb0); % TODO: what should sbmin, sbmax, alpha, and sb0 be 
%     sb = log(1 + exp(alpha * (sb - 0.01)))/alpha;
    % Artificial Bulk viscosity
    avb = (kb.*hm./(porder)) .* (abs(uv) + abs(vv) + abs(c_star)) .* sb;
    
    % Assign artificial viscosities
    avField(1) = avb;  %  bulk
end

function m = mass(u, q, w, v, x, t, mu, eta)
    nspecies = 5;
    ndim = 2;
    m = sym(ones(nspecies + ndim + 1, 1));
end

function f = flux(u, q, w, v, x, t, mu, eta)
    nspecies = 5;
    nd = 2;
    % % Clipper utilities
    % alpha = 1.0e3;
    % rmin = 0.0;

    fi = fluxinviscid(u, q, w, v, x, t, mu, eta, nspecies, nd);
    % fv = fluxviscous(u, q, w, v, x, t, mu, eta)
    fv = fluxlaplacian(u, q, w, v, x, t, mu, eta, nspecies, nd);
    f = fi - fv;
end

function fv = fluxlaplacian(u, q, w, v, x, t, mu, eta, nspecies, ndim)
    rmin = 0.0;
    alpha = 1.0e12;
    eps_av = v(1);
%     hm = v(2);
    %TODO: CHECK GAS OUTPUTS TO INPUTS USED BY HTR GROUP 
    %TODO: CHECK NONDIM USED BY HTR GROUP !
    r_vec = u(1:nspecies);
    dr = dlmax(r_vec-rmin,alpha); % Regularize derivative 
    fv(1:nspecies) = eps_av*-q(1:nspecies) .* dr; % drho_i
    fv(nspecies+1) = eps_av*-q(nspecies+1); % drhou
    fv(nspecies+2) = eps_av*-q(nspecies+2); % drhov
    fv(nspecies+3) = eps_av*-q(nspecies+3); % drhoE
end

function fi = fluxinviscid(u, q, w, v, x, t, mu, eta, nspecies, ndim)
    nenergy = 1;
    % ndim = 2;
    fi = sym(zeros(nspecies+ndim+1),2);

    rho_i = sym(zeros(nspecies,1));
    rho = sym(0);

    rmin = 0.0;
    alpha = 1.0e8;

    % Conservative Variables
    for ispecies = 1:nspecies
        rho_i(ispecies) = rmin + lmax(u(ispecies)-rmin,alpha); %subspecies density
        rho = rho + rho_i(ispecies); %total mixture density
    end

    rhou = u(nspecies+1);
    if ndim == 2
        rhov = u(nspecies+2)
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



function fv = fluxviscous(u, q, w, v, x, t, mu, eta, nspecies, ndim)
    %TODO for thursday:
    %      - double check with papers: consistent with may, gnoffo. Maybe typo in su2mpp? 
    %      - double check with NS code (TODO: REALLY CONFUSED ABOUT SIGN OF STRESS TENSOR)
    %      - add AV options: tomorrow, be clear about shear, bulk, thermal, and how density should come into this (let's discuss)
    %      - try to generate 
    rmin = 0.0;
    alpha = 1e12;
    nenergy=1;
    ncu = nspecies + ndim + nenergy;

    Pr = mu(18);
    Re = mu(19);

    % Data layout
    %   [D_1, ..., D_ns, h_1, ..., h_ns, mu, kappa, dTdri, dTdrhou, dTdrhoE]
    p = eta(1);
    D_vec = eta(2:nspecies+1); 
    h_vec = eta(nspecies+2:2*nspecies+1);
    viscshear = eta(2*nspecies+2);
    kappa = eta(2*nspecies+3);
    dT_drho_i = eta(2*nspecies+4:3*nspecies+3);
    dT_drhoe = eta(3*nspecies+4);

    rho_i = rmin + lmax(u(1:nspecies) - rmin, alpha);

    rhou = u(nspecies+1);
    rhov = u(nspecies+2);
    rhoE = u(nspecies+ndim+1);
    
    drho_dx_i = -q(1:nspecies) .* dlmax(u(1:nspecies)-rmin,alpha)
    drhou_dx  = -q(nspecies+1);
    drhov_dx  = -q(nspecies+2)
    drhoE_dx  = -q(nspecies+ndim+1);
    drho_dy_i = -q(ncu+1:ncu+nspecies) .* dlmax(q(ncu+1:ncu+nspecies)-rmin,alpha)
    drhou_dy  = -q(ncu+nspecies+1);
    drhov_dy  = -q(ncu+nspecies+2)
    drhoE_dy  = -q(ncu+nspecies+ndim+1);

    % Some useful derived quantities 
    rho = sum(rho_i);
    drho_dx = sum(drho_dx_i);
    drho_dy = sum(drho_dy_i);
    rho_inv = 1.0 / rho;
    drho_dx_inv = 1.0 / drho_dx;
    drho_dy_inv = 1.0 / drho_dy
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
    kappa = kappa + avk; % thermal conductivity
    mu_d = mu_d + avs;   % dynamic viscosity
    beta = 0 + avb;  % bulk viscosity 
    D_vec = D_vec + avDvec; % Species diffusion velocities

    %%%%%%%% Calculation of J_i
    Y_i = rho_i ./ rho;
    dY_dx_i = (drho_dx_i * rho - rho_i * drho_dx) * rho_inv * rho_inv;
    dY_dy_i = (drho_dy_i * rho - rho_i * drho_dy) * rho_inv * rho_inv;
 
    J_i_x = -rho * D_vec .* dY_dx_i + rho_i .* sum(D_vec .* dY_dx_i);
    J_i_y = -rho * D_vec .* dY_dy_i + rho_i .* sum(D_vec .* dY_dy_i);

    %%%%%%%% Stress tensor tau
    % txx = viscshear * 4.0/3.0 * du_dx + viscbulk * du_dx; %TODO: check sign
    txx = mu_d * 2.0/3.0 * (2 * du_dx - dv_dy) / Re + beta * (du_dx + dv_dy);
    txy = mu_d * (du_dy + dv_dx) / Re;
    tyy = mu_d * 2.0/3.0 * (2 * dv_dy - du_dx) / Re + beta * (du_dx + dv+dy)
    %%%%%%%% Final viscous fluxes
    for ispecies = 1:nspecies
        fv(ispecies,1) = -J_i_x(ispecies);
        fv(ispecies,2) = -J_i_y(ispecies);
    end

    fv(nspecies + 1, 1) = txx;
    fv(nspecies + 2, 1) = txy;
    fv(nspecies + ndim + 1,1) = uv * txx + vv * txy - (sum(h_vec.*J_i_x) - kappa * dT_dx / (Re * Pr));
    fv(nspecies + 1, 2) = txy;
    fv(nspecies + 2, 2) = tyy
    fv(nspecies + ndim + 1,2) = uv * txy + vv * tyy - (sum(h_vec.*J_i_y) - kappa * dT_dy / (Re * Pr));
end


function s = source(u, q, w, v, x, t, mu, eta)
    
    nspecies = 5;
    rmin = 0.0;
    alpha = 1.0e12;
    rvec = u(1:nspecies);
    s(1:nspecies) = rmin + lmax(rvec-rmin,alpha);
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    nspecies = 5;
    ndim = 2;
    ub = sym(zeros(nspecies+2, 2));
%%%%%% C-CODE MANUALLY WRITTEN
%     ub(:,1) = u(:);
%     ub(:,2) = u(:); 

    uinflow = mu(1:nspecies+2+1);

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
    fb = sym(zeros(nspecies+2, 2));

    f = flux(u, q, w, v, x, t, mu, eta);
    fn = f(:,1)*n(1) + f(:,2)*n(2) + tau.*(u-uhat);

    fadiabatic = fi;
    fadiabatic(1:nspecies) = 0.0;
    fadiabatic(nspecies+ndim+2) 0.0;
    fb(:,1) = fn;
    fb(:,2) = fn;
    fb(:,3) = fadiabatic;
end

function u0 = initu(x, mu, eta)
    nspecies = 5;
    u0 = mu(1:nspecies+2+1);
end

function o = output(u, q, w, v, x, t, mu, eta)
    
    % nspecies = 5;

    % o = sym(zeros(10,1));
    % nspecies = 5;
    % alpha = 1.0e3;
    % rmin = 0.0;

    % f = sym(zeros(nspecies + 2,1));
    % rho_i = sym(zeros(5,1));
    % rho = sym(0);

    % % Conservative Variables
    % for ispecies = 1:nspecies
    %     rho_i(ispecies) = rmin + lmax(u(ispecies)-rmin,alpha); %subspecies density
    %     rho = rho + rho_i(ispecies); %total mixture density
    % end

    % rhou = u(nspecies+1);
    % rhoE = u(nspecies+2);

    % rhoinv = 1.0 / rho;
    % uv = rhou * rhoinv; %velocity
    % E = rhoE*rhoinv; %energy

    % o(1) = u(1);
    % o(2) = u(2);
    % o(3) = u(3);
    % o(4) = u(4);
    % o(5) = u(5);
    o = sensor_outputs(u, q, w, v, x, t, mu, eta);
end

% function dout = dlmax(x, alpha)
%     dout = atan(alpha*(x))/pi + (alpha*(x))./(pi*(alpha^2*(x).^2 + 1)) + 1/2;
% end