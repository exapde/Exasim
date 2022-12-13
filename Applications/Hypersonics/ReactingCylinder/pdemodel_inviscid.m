 
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
    f = avfieldLaplacian2d(u, q, w, v, x, t, mu, eta, nspecies, nd);
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
    fl = fluxlaplacian(u, q, w, v, x, t, mu, eta, nspecies, nd);
    f = fi - fl;
end

function fv = fluxlaplacian(u, q, w, v, x, t, mu, eta, nspecies, ndim)
    rmin = 0.0;
    alpha = 1.0e12;
    eps_av = v(1);
    nenergy = 1;
    ncu = nspecies + ndim + nenergy;
%     hm = v(2);
    %TODO: CHECK GAS OUTPUTS TO INPUTS USED BY HTR GROUP 
    %TODO: CHECK NONDIM USED BY HTR GROUP !
    r_vec = u(1:nspecies);
    dr = dlmax(r_vec-rmin,alpha); % Regularize derivative 

    drho_dx_i = -q(1:nspecies).*dr; %.* dlmax(u(1:nspecies)-rmin,alpha);
    drhou_dx  = -q(nspecies+1);
    drhov_dx  = -q(nspecies+2);
    drhoE_dx  = -q(nspecies+ndim+1);
    drho_dy_i = -q(ncu+1:ncu+nspecies) .*dr;% .* dlmax(q(ncu+1:ncu+nspecies)-rmin,alpha);
    drhou_dy  = -q(ncu+nspecies+1);
    drhov_dy  = -q(ncu+nspecies+2);
    drhoE_dy  = -q(ncu+nspecies+ndim+1);

    fv(1:nspecies,1) = eps_av*drho_dx_i; % dxrho_i
    fv(nspecies+1,1) = eps_av*drhou_dx; % dxrhou
    fv(nspecies+2,1) = eps_av*drhov_dx; % dxrhov
    fv(nspecies+3,1) = eps_av*drhoE_dx; % dxrhoE
    fv(1:nspecies,2) = eps_av*drho_dy_i; % dxrho_i
    fv(nspecies+1,2) = eps_av*drhou_dy; % dxrhou
    fv(nspecies+2,2) = eps_av*drhov_dy; % dxrhov
    fv(nspecies+3,2) = eps_av*drhoE_dy; % dxrhoE
end

function fi = fluxinviscid(u, q, w, v, x, t, mu, eta, nspecies, ndim)
    nenergy = 1;
    % ndim = 2;
    fi = sym(zeros(nspecies+ndim+1,2));

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

function s = source(u, q, w, v, x, t, mu, eta)
    
    nspecies = 5;
    rmin = 0.0;
    alpha = 1.0e12;
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
    nx = n(1);
    ny = n(2);
% % 
    uadiabatic(nspecies+1) = u(nspecies+1) - nx * (u(nspecies+1)*nx + u(nspecies+2)*ny);
    uadiabatic(nspecies+2) = u(nspecies+2) - ny * (u(nspecies+1)*nx + u(nspecies+2)*ny);

    uadiabatic = 0.5*(uadiabatic + u);
    ub(:,1) = uinflow;
    ub(:,2) = uoutflow;
    ub(:,3) = uadiabatic;
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    nspecies = 5;
    ndim = 2;
    fb = sym(zeros(nspecies+ndim+1, 3));

    f = flux(u, q, w, v, x, t, mu, eta); % maybe should be uhat here? 
    fn = f(:,1)*n(1) + f(:,2)*n(2) + tau.*(u-uhat);

%     OH maybe I need to change 
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

    o(1:8) = u(1:8);
    o(1) = v(1);
%     o = sensor_outputs(u, q, w, v, x, t, mu, eta);
end

% function dout = dlmax(x, alpha)
%     dout = atan(alpha*(x))/pi + (alpha*(x))./(pi*(alpha^2*(x).^2 + 1)) + 1/2;
% end