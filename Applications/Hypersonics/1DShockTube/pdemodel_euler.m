 
function pde = pdemodel_inflow
    pde.mass = @mass;
    pde.flux = @flux;
    pde.source = @source;
    pde.fbou = @fbou;
    pde.ubou = @ubou;
    pde.initu = @initu;
end

function m = mass(u, q, w, v, x, t, mu, eta)
    nspecies = 1;
    m = sym(ones(nspecies + 2, 1));
%     m = sym([1.0; 1.0; 1.0; 1.0; 1.0]); 
end

function f = flux(u, q, w, v, x, t, mu, eta)

    nspecies = 1;

    f = sym(zeros(nspecies + 2,1));
    
    rho = u(1);
    rhou = u(2);
    rhoE = u(3);

    rhoinv = 1.0 / rho;
    uv = rhou * rhoinv; %velocity
    E = rhoE * rhoinv; %energy

    gamma = 1.4;
    gammam1 = gamma-1;
    p = gammam1 * (rhoE - 0.5 * rho * uv^2);
    H = E+p*rhoinv; %enthalpy

    f(1) = rhou;
    f(nspecies + 1) = rhou * uv + p;
    f(nspecies + 2) = rhou * H;
end

function s = source(u, q, w, v, x, t, mu, eta)
    nspecies = 1;
    s = sym(zeros(nspecies + 2,1));
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    nspecies = 1;
    ub = sym(zeros(nspecies+2, 2));

%     rho = abs(u(1));
%     rhou = u(2);
%     rhoE = u(3);
% 
%     rhoinv = 1.0 / rho;
%     uv = rhou * rhoinv; %velocity
%     E = rhoE * rhoinv; %energy
% 
%     gamma = 1.4;
%     gammam1 = gamma-1;
%     p = gammam1 * (rhoE - 0.5 * rho * uv^2);
%     H = E+p*rhoinv; %enthalpy
%     a = sqrt(gamma * p / rho);
%     un = uv*n; %variable in case we u * n; dont think it's necessary for 1d
% 
%     % From Toro p. 90
%     K = [1.0,            1.0,       1.0;
%          un-a,           un,        un+a;
%          H - un * a,  0.5 * un^2,  H + un*a];
% 
%     % From Toro p. 276
%     Kinv = (gammam1)/(2*a^2) * [0.5 * un^2 + (un * a)/(gammam1), -un - a/gammam1, 1.0;...
%                               2*a^2/gammam1 - un^2,                   2*un,     -2.0;...
%                               0.5*un^2 - un*a/gammam1,       a/gammam1 - un,    1.0];
%     
%     % Taken from cylinder pdemodel; picks sign of eigenvalue
%     Lambda = [tanh(1e3 * (un - a)), 0,        0;
%                 0,              tanh(1e3*un), 0; 
%                 0,               0,         tanh(1e3*(un+a))];
%     A = K * Lambda * Kinv;
% 
%     Minf = eta(4);
%     uinf = [1.0; 1.0; 1.0/(gamma * gammam1 * Minf^2) + 0.5];

    uin = uinflow(u, mu, eta);
%     uin = 0.5 * ((u + uinf) + A * (u - uinf));
    ub(:,1) = uin;

    uout = uoutflow(u, mu, eta);
    ub(:,2) = uout;
end  

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    nspecies = 1;
    fb = sym(zeros(nspecies+2, 2));

%     uin = uinflow(u, mu, eta);
%     fin = flux(uin, q, w, v, x, t, mu, eta);
    f = flux(u, q, w, v, x, t, mu, eta);
    fb(:,1) = f *n(1) + tau(1)*(u-uhat); 
    fb(:,2) = f *n(1) + tau(1)*(u-uhat); 
%     uout = uoutflow(u, mu, eta);
%     fout = flux(uout, q, w, v, x, t, mu, eta);
%     fb(:,2) = fout*n(1) + tau(1)*(uout-uhat);
end

function u0 = initu(x, mu, eta)
    nspecies = 1;
%     u0 = uinflow(mu,mu,eta);
%     u0 = uoutflow(mu,mu,eta);
%     u0 = mu(1:3);
    gamma = 1.4;
    gammam1 = gamma - 1.0;
    Minf = eta(4);
    u0 = [1.0; 1.0; 1.0/(gamma * gammam1 * Minf^2) + 0.5];
end

function uin = uinflow(u, mu, eta)
    nspecies = 1;
    % Typical subsonic inflow: given u, rho, grab E

    uin = [1; 1; u(3)];
    % Equivalent to nonequil formulation; given u, T, X=1
%     rho = u(1);
%     rhou = u(2);
%     rhoE = u(3);
% 
%     rhoinv = 1.0/rho;
%     uv = rhou * rhoinv;
%     gamma = 1.4;
%     gammam1 = gamma - 1.0;
% 
%     p = gammam1 * (rhoE - 0.5 * rho * uv^2);
% 
%     T_in = eta(5);
%     u_in = eta(2);
% 
%     rho_in = p / (gammam1 * T_in);
%     rhou_in = rho_in * u_in;
%     rhoE_in = p/gammam1 + 0.5 * rho_in * u_in^2;
% 
%     uin = [rho_in; rhou_in; rhoE_in];

%     rhou = eta(2)*rho;
% %     rhou = eta(2)*rho;
%     rhoE = u(3);
% 
%     rhoinv = 1.0 / rho;
%     uv = rhou * rhoinv; %velocity
%     E = rhoE * rhoinv; %energy
% 
%     gamma = 1.4;
%     gammam1 = gamma-1;
%     p = gammam1 * (rhoE - 0.5 * rho * uv^2);
% 
%     uv_inflow = mu(2);
%     uin = [rho; rhou; rhoE];
end

function uout = uoutflow(u, mu, eta)
    nspecies = 1;

    gamma = 1.4;
    gammam1 = gamma-1;
    p = eta(3);

    rho = u(1);
    rhou = u(2);
    rhoinv = 1/rho;
    uv = rhou * rhoinv;

    rhoE = p/gammam1 + 0.5 * rho * uv^2;
    uout = [rho; rhou; rhoE];
end
