 
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
    nd = 1;
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
    % rv = udg(3);
    rE = u(nspecies+nd+1);

    rx = sum(-q(1:nspecies));
    rux = -q(nspecies+1);

    % Regularization 
    rvec = rmin + lmax(rvec-rmin,alpha); % need to double check regularizatino here
    
    % Simple derived working quantities
    r = sum(rvec);
    r1 = 1./r;
    uv = ru.*r1;
    E = rE * r1; %energy
    H = E + p*r1; %enthalpy
    H = Hmin + lmax(H - Hmin, alpha);

    % Critical speed of Sound
    c_star = sqrt((2.*gam1.*H) ./ (gam+1)); %TODO: this is fine for 1 temp but not 2

    % Computing derivatives for the sensors
    ux = (rux - rx.*uv).*r1;
    div_v = (ux); %TODO: check signs. Actually this probably is important, we want to make sure it's applied in negative v. pos dilitation
                  %      pretty sure my signs are okay. There is something strange about the code I was given. Truly don't understand the signs 
    % limit  divergence and vorticity
    sigm = 1e4;
    div_v = limiting(div_v,-sigm,sigm,alpha,-sigm);

    % % Dilatation Sensor sb
    sb = - (hm./porder) .* (div_v./c_star);
    sb = limiting(sb,sbmin,sbmax,alpha,sb0); % TODO: what should sbmin, sbmax, alpha, and sb0 be 
%     sb = log(1 + exp(alpha * (sb - 0.01)))/alpha;
    % Artificial Bulk viscosity
    avb = (kb.*hm./(porder)) .* (abs(uv) + abs(c_star)) .* sb;
    
    % Assign artificial viscosities
    avField(1) = avb;  %  bulk
end

function m = mass(u, q, w, v, x, t, mu, eta)
    nspecies = 5;
    m = sym(ones(nspecies + 2, 1));
end

function f = flux(u, q, w, v, x, t, mu, eta)
    nspecies = 5;
    nd = 1;
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
    fv(nspecies+2) = eps_av*-q(nspecies+2); % drhoE
end

function fi = fluxinviscid(u, q, w, v, x, t, mu, eta, nspecies, ndim)
    nenergy = 1;

    rho_i = sym(zeros(nspecies,1));
    rho = sym(0);

    rmin = 0.0;
    alpha = 1.0e12;

    % Conservative Variables
    for ispecies = 1:nspecies
        rho_i(ispecies) = rmin + lmax(u(ispecies)-rmin,alpha); %subspecies density
        rho = rho + rho_i(ispecies); %total mixture density
    end

    rhou = u(nspecies+1);
    % if ndim == 2
    %     rhov = u(nspecies+2)
    % end
    rhoE = u(nspecies+ndim+1);
    if nenergy == 2
        rhoev = u(nspecies+ndim+nenergy);
    end

    rhoinv = 1.0 / rho;
    uv = rhou * rhoinv; %velocity
    E = rhoE * rhoinv; %energy

    p = eta(1); % MANUALLY MODIFIED
    H = E + p*rhoinv; %enthalpy

    % Fluxes
    for ispecies = 1:nspecies
        fi(ispecies) = rho_i(ispecies) * uv;
    end
    fi(nspecies + 1) = rhou * uv + p;
    fi(nspecies + 2) = rhou * H;
end



function fv = fluxviscous(u, q, w, v, x, t, mu, eta)
    %TODO for thursday:
    %      - double check with papers: consistent with may, gnoffo. Maybe typo in su2mpp? 
    %      - double check with NS code (TODO: REALLY CONFUSED ABOUT SIGN OF STRESS TENSOR)
    %      - add AV options: tomorrow, be clear about shear, bulk, thermal, and how density should come into this (let's discuss)
    %      - try to generate 
    nspecies = 5;
    nenergy = 1;
    ndim = 1;
    alpha = 1e9;

    Pr = mu(18);
    Re = mu(19);

    avb = v(1);
    avs = v(2);
    avk = v(3);
    avDvec = v(4:3+nspecies);

    % Data layout
    %   [D_1, ..., D_ns, h_1, ..., h_ns, mu, kappa, dTdri, dTdrhou, dTdrhoE]
    D_vec = eta(1:nspecies); 
    h_vec = eta(nspecies+1:2*nspecies);
    viscshear = eta(2*nspecies+1);
    kappa = eta(2*nspecies+2);
    dT_drho_i = eta(2*nspecies+3:3*nspecies+2);
    dT_drhoe = eta(3*nspecies+3);

    rho_i = rmin + lmax(u(1:nspecies) - rmin, alpha);
    drho_dx_i = -q(1:nspecies) .* dlmax(u(1:nspecies)-rmin,alpha)
    rho = sum(rho_i);
    drho_dx = sum(drho_dx_i);
    rhou     = u(nspecies+1);
    drhou_dx = -q(nspecies+1);
    rhoE     = u(nspecies+ndim+1);
    drhoE_dx = -q(nspecies+ndim+1);

    % Some useful derived quantities 
    rho_inv = 1.0 / rho;
    drho_dx_inv = 1.0 / drho_dx;
    uv = rhou * rho_inv; %velocity
    dudx = (drhou_dx - drho_dx*uv)*rho_inv;
    % re = rhoE - rho * (0.5 * uTu); 
    uTu2 = 0.5 * uv * uv;
    duTu2_dx = uv; 
    dre_drho = -uTu2;
    dre_duTu2 = -rho;
    dre_drhoE = 1.0;
    dre_dx = dre_drho * drho_dx + dre_duTu2 * duTu2_dx + dre_drhoE * drhoE_dx;
    Tx = sum(dT_drho_i .* drho_dx_i) +  dT_drhoe * dre_dx;
    
    %%%%%%%% Modify appropriate coefficients 
    kappa = kappa + avk;
    viscshear = viscshear + avs;
    viscbulk = 0 + avb;
    D_vec = D_vec + avDvec;

    %%%%%%%% Calculation of J_i
    % for ispecies = 1:nspecies
    %     Y_i(ispecies) = rho_i(ispecies) / rho;
    %     dYdx_i(ispecies) = (q(ispecies) * rho - rho_i(ispecies) * drhodx) * drhodxinv^2;
    % end
    Y_i = rho_i ./ rho;
    dY_dx_i = (drho_dx_i * rho - rho_i * drho_dx) * rho_inv * rho_inv;
    
    % for ispecies = 1:nspecies
        % J_i(ispecies) = -rho * D_vec(ispecies) * dYdx_i(ispecies) + rho_i(ispecies) * sum(D_vec .* dYdx_i);
    % end
    J_i = -rho * D_vec .* dY_dx_i + rho_i .* sum(D_vec .* dY_dx_i);

    %%%%%%%% Stress tensor tau
    txx = 1.0 / Re * viscshear * 4.0/3.0 * du_dx + viscbulk * du_dx; %TODO: check sign

    %%%%%%%% Final viscous fluxes
    for ispecies = 1:nspecies
        f(ispecies) = -J_i(ispecies);
    end

    fv(nspecies + 1) = txx;
    fv(nspecies + 2) = uv * txx - (sum(h_vec.*J_i) - kappa * dT_dx * 1.0 / (Re * Pr));
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
    ub = sym(zeros(nspecies+2, 2));
%%%%%% C-CODE MANUALLY WRITTEN
%     ub(:,1) = u(:);
%     ub(:,2) = u(:); 
    ub(:,1) = initu(x, mu, eta);
    ub(:,2) = initu(x, mu, eta);
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    nspecies = 5;
    fb = sym(zeros(nspecies+2, 2));

    f = flux(u, q, w, v, x, t, mu, eta);
    fn = f*n(1);
%%%%% C-CODE MANUALLY MODIFIED
    for i = 1:(nspecies+2)
        fb(i,1) = fn(i) + tau(i)*(u(i)-uhat(i));
    end
    fb(:,2) = fb(:,1);
end

function u0 = initu(x, mu, eta)
    nspecies = 5;

    % pde.physicsparam = [rho_post(:)', rhou_post, rhoE_post, rho_equil(:)', rhou_equil, rhoE_equil,   Xi_inflow,    u_inflow, T_inflow, p_outflow];
    %                    %  1:N          N+1         N+2       N+3:2*N+2     2*N+3       2*N+4       2*N+5:3*N+4     3*N+5     3*N+6      3*N+7
    % pde.externalparam = [rho_scale, u_scale, rhoe_scale]; 
    % oh let's just modify this manually 
    rho_L_vec = mu(1:nspecies);
    rhou_L = mu(nspecies + 1);
    rhoE_L = mu(nspecies + 2);

    rho_R_vec = mu(nspecies+3:2*nspecies+2);
    rhou_R = mu(2*nspecies + 3);
    rhoE_R = mu(2*nspecies + 4); 
    %%% TODO: Cinnella scales this with right region. And uses speed of sound. BUT WHY IS ENERGY IN THE RIGHT REGION NEGATIVE??? 
    rho_scale = eta(1);
    u_scale = eta(2);
    rhoe_scale = eta(3);

    rho_R_vec_nondim = rho_R_vec/rho_scale;
    rhou_R_nondim = rhou_R/(u_scale * rho_scale);
    rhoE_R_nondim = rhoE_R/rhoe_scale;

    rho_L_vec_nondim = rho_L_vec/rho_scale;
    rhou_L_nondim = rhou_L/(u_scale * rho_scale);
    rhoE_L_nondim = rhoE_L/rhoe_scale;

%     u0 = [rho_R_vec_nondim(:) + rho_L_vec_nondim(:); rhou_R_nondim + rhou_L_nondim; rhoE_R_nondim + rhoE_L_nondim];
    stepsize = 0.001; 
    midpoint = 0.5;
    u0(1) = smoothstep_down(x-midpoint, rho_L_vec_nondim(1), rho_R_vec_nondim(1), stepsize);
    u0(2) = smoothstep_down(x-midpoint, rho_L_vec_nondim(2), rho_R_vec_nondim(2), stepsize);
    u0(3) = smoothstep_down(x-midpoint, rho_L_vec_nondim(3), rho_R_vec_nondim(3), stepsize);
    u0(4) = smoothstep_up(x-midpoint, rho_L_vec_nondim(4), rho_R_vec_nondim(4), stepsize);
    u0(5) = smoothstep_up(x-midpoint, rho_L_vec_nondim(5), rho_R_vec_nondim(5), stepsize);
    u0(nspecies+1) = 0.0;
    u0(nspecies+2) = smoothstep_down(x-midpoint, rhoE_L_nondim, rhoE_R_nondim, stepsize);

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
    hm = v(2);
    o = sensor_outputs(u, q, w, v, x, t, mu, eta, hm);
end

% function dout = dlmax(x, alpha)
%     dout = atan(alpha*(x))/pi + (alpha*(x))./(pi*(alpha^2*(x).^2 + 1)) + 1/2;
% end