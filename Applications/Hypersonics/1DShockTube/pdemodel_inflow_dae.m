 
function pde = pdemodel_inflow
    pde.mass = @mass;
    pde.flux = @flux;
    pde.source = @source;
    pde.fbou = @fbou;
    pde.ubou = @ubou;
    pde.initu = @initu;
    pde.sourcew = @sourcew;
    pde.initw = @initw;
end

function m = mass(u, q, w, v, x, t, mu, eta)
    nspecies = 5;
    m = sym(ones(nspecies + 2, 1));
%     m = sym([1.0; 1.0; 1.0; 1.0; 1.0]); 
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

    % % Mutation outputs
	% rho_inf = eta(1);
    % u_inf = eta(2);

    p = w(nspecies + 1); %pressure, assume nondimensionalized

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
        s(ispecies) = w(ispecies);
    end
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    nspecies = 5;
    ub = sym(zeros(nspecies+2, 2));
%     [rho_pre, rhou_pre, rhoE_pre, rho_post, rhou_postshock, rhoE_post]
%     rhoL = mu(1:nspecies);
%     rhouL = mu(nspecies+1);
%     rhoEL = mu(nspecies+2);
% 
%     rho_post = mu(nspecies+3:2*nspecies+2);
%     rhou_post = mu(2*nspecies+3);
%     rhoE_post = mu(2*nspecies+4);

%     ub(:,1) = [rhoL(:); rhouL; rhoEL];
    
%     ub(:,1) = mu(8:14);
%     ub(:,1) = mu(1:nspecies+2);
%     ub(:,1) = u(:);
%     p_outflow = mu(nspecies+3);
%     ub(:,2) = [u(1:nspecies);u(nspecies+1);p_outflow];
%%%%%% C-CODE MANUALLY WRITTEN
    ub(:,1) = u(:);
    ub(:,2) = u(:); 
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    nspecies = 5;
    fb = sym(zeros(nspecies+2, 2));

    f = flux(u, q, w, v, x, t, mu, eta);

%%%%% C-CODE MANUALLY WRITTEN
    fb(:,1) = f*n(1) + tau(1)*(u-uhat);
    fb(:,2) = f*n(1) + tau(1)*(u-uhat);
end

function u0 = initu(x, mu, eta)
    % To be manually modified
    nspecies = 5;
    % TODO:
%     u0 = sym(zeros(nspecies+2,1));
% 
%     rho_pre = mu(1:nspecies);
%     rhou_pre = mu(nspecies+1);
%     rhoE_pre = mu(nspecies+2);
% 
%     rho_post = mu(nspecies+3:2*nspecies+2);
%     rhou_post = mu(2*nspecies+3);
%     rhoE_post = mu(2*nspecies+4);
%     u0 = mu(nspecies + 3:2*nspecies + 2);
    % u0 = mu(1:nspecies+2);

%     pde.physicsparam = [rho_post(:)', rhou_post, rhoE_post, rho_equil(:)', rhou_equil, rhoE_equil,   Xi_inflow,    u_inflow, T_inflow, p_outflow]
%                           1:N          N+1         N+2       N+3:2*N+2     2*N+3       2*N+4       2*N+5:3*N+4     3*N+5     3*N+6      3*N+7

    rho_post_vec = mu(1:nspecies);
    rhou_post = mu(nspecies + 1);
    rhoE_post = mu(nspecies + 2);

    rho_inf = eta(1);
    u_inf = eta(2);
    rhoe_inf = eta(3);

    rho_equil_vec = mu(nspecies+3:2*nspecies+2);
    rhou_equil = mu(2*nspecies + 3);
    rhoE_equil = mu(2*nspecies + 4);

    rho_equil_vec_nondim = rho_equil_vec/rho_inf;
    rhou_equil_nondim = rhou_equil/(u_inf * rho_inf);
    rhoE_equil_nondim = rhoE_equil/rhoe_inf;

    rho_post_vec_nondim = rho_post_vec/rho_inf;
    rhou_post_nondim = rhou_post/(u_inf * rho_inf);
    rhoE_post_nondim = rhoE_post/rhoe_inf;

%     u0 = [rho_equil_vec_nondim(:); rhou_equil_nondim; rhoE_equil_nondim];
%     u0 = [rho_post_vec_nondim(:); rhou_post_nondim; rhoE_post_nondim];
%     for i = 1:nspecies
%         u0(i) = smoothstep(x, rho_post(i), rho_pre(i), 1/1000);
%     end
    stepsize = 0.01;
    u0(1) = smoothstep_up(x-stepsize, rho_post_vec(1)/rho_inf, rho_equil_vec_nondim(1), stepsize);
    u0(2) = smoothstep_up(x-stepsize, rho_post_vec(2)/rho_inf, rho_equil_vec_nondim(2), stepsize);
    u0(3) = smoothstep_up(x-stepsize, rho_post_vec(3)/rho_inf, rho_equil_vec_nondim(3), stepsize);
    u0(4) = smoothstep_down(x-stepsize, rho_post_vec(4)/rho_inf, rho_equil_vec_nondim(4), stepsize);
    u0(5) = smoothstep_down(x-stepsize, rho_post_vec(5)/rho_inf, rho_equil_vec_nondim(5), stepsize);
    u0(nspecies+1) = smoothstep_down(x-stepsize, rhou_post/(rho_inf*u_inf), rhou_equil_nondim, stepsize);
    u0(nspecies+2) = smoothstep_down(x-stepsize, rhoE_post/rhoe_inf, rhoE_equil_nondim, stepsize);
%     u0 = [rhoL(:) + rhoR(:); rhouL + rhouR; rhoEL + rhoER];
end

function sw = sourcew(u, q, w, v, x, t, mu, eta)
    % MODIFIED MANUALLY
    nspecies = 5;

    sw = sym(zeros(nspecies + 2,1));
    rho_i = sym(zeros(5,1));
    rho = sym(0);

    % Conservative Variables
    for ispecies = 1:nspecies
        rho_i(ispecies) = u(ispecies); %subspecies density
        rho = rho + rho_i(ispecies); %total mixture density
    end

    rhou = u(nspecies+1);
    rhoE = u(nspecies+2);

    rhoinv = 1.0 / rho;
    u = rhou * rhoinv; %velocity
    E = rhoE*rhoinv; %energy

    % Mutation inputs
    rhoe = rhoE - 0.5 * rho * u^2;
    for i = 1:nspecies
        sw(i) = rho_i(i);
    end
    sw(nspecies + 1) = rhoe;
end

function w0 = initw(x, mu, eta)
    % MODIFIED MANUALLY: could be specified but it's a nice test of the initial conditions. 
    nspecies = 5;
    w0 = sym(zeros(nspecies + 1,1));
end