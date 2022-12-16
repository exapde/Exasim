 
function pde = pdemodel_inflow
    pde.mass = @mass;
    pde.flux = @flux;
    pde.source = @source;
    pde.fbou = @fbou;
    pde.ubou = @ubou;
    pde.initu = @initu;
    pde.output = @output;
end

function m = mass(u, q, w, v, x, t, mu, eta)
    nspecies = 5;
    m = sym(ones(nspecies + 2, 1));
end

function f = flux(u, q, w, v, x, t, mu, eta)

    nspecies = 5;
    alpha = 1.0e3;
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

    p = eta(4); % MANUALLY MODIFIED
    H = E+p*rhoinv; %enthalpy

    % Fluxes
    for ispecies = 1:nspecies
        f(ispecies) = rho_i(ispecies) * uv;
    end
    f(nspecies + 1) = rhou * uv + p;
    f(nspecies + 2) = rhou * H;
%     f(nspecies + 2) = (rhoE + p)*u;
end

function s = source(u, q, w, v, x, t, mu, eta)
    
    nspecies = 5;

    s = sym(zeros(nspecies + 2,1));
    omega_i = sym(zeros(nspecies,1));

    for ispecies = 1:nspecies
        s(ispecies) = 1; % MANUALLY MODIFIED
    end
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    nspecies = 5;
    ub = sym(zeros(nspecies+2, 2));
%%%%%% C-CODE MANUALLY WRITTEN
    ub(:,1) = u(:);
    ub(:,2) = u(:); 
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
    
    rho_post_vec = mu(1:nspecies);
    rhou_post = mu(nspecies + 1);
    rhoE_post = mu(nspecies + 2);

    rho_equil_vec = mu(nspecies+3:2*nspecies+2);
    rhou_equil = mu(2*nspecies + 3);
    rhoE_equil = mu(2*nspecies + 4);

    rho_scale = eta(1);
    u_scale = eta(2);
    rhoe_scale = eta(3);

    rho_equil_vec_nondim = rho_equil_vec/rho_scale;
    rhou_equil_nondim = rhou_equil/(u_scale * rho_scale);
    rhoE_equil_nondim = rhoE_equil/rhoe_scale;

    rho_post_vec_nondim = rho_post_vec/rho_scale;
    rhou_post_nondim = rhou_post/(u_scale * rho_scale);
    rhoE_post_nondim = rhoE_post/rhoe_scale;
% Smooth step
%     stepsize = 0.01; 
%     u0(1) = smoothstep_up(x-stepsize, rho_post_vec_nondim(1), rho_equil_vec_nondim(1), stepsize);
%     u0(2) = smoothstep_up(x-stepsize, rho_post_vec_nondim(2), rho_equil_vec_nondim(2), stepsize);
%     u0(3) = smoothstep_up(x-stepsize, rho_post_vec_nondim(3), rho_equil_vec_nondim(3), stepsize);
%     u0(4) = smoothstep_down(x-stepsize, rho_post_vec_nondim(4), rho_equil_vec_nondim(4), stepsize);
%     u0(5) = smoothstep_down(x-stepsize, rho_post_vec_nondim(5), rho_equil_vec_nondim(5), stepsize);
%     u0(nspecies+1) = smoothstep_down(x-stepsize, rhou_post_nondim, rhou_equil_nondim, stepsize);
%     u0(nspecies+2) = smoothstep_up(x-stepsize, rhoE_post_nondim, rhoE_equil_nondim, stepsize);
    u0 = [rho_equil_vec_nondim(:); rhou_equil_nondim; rhoE_equil_nondim];
end

function o = output(u, q, w, v, x, t, mu, eta)
    
    nspecies = 5;

    o = sym(zeros(2,1));
    nspecies = 5;
    alpha = 1.0e3;
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

    o(1) = u(1);
    o(2) = u(2);
end