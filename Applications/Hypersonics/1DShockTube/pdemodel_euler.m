 
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

    uin = uinflow(u, mu, eta);
    ub(:,1) = uin;

    uout = uoutflow(u, mu, eta);
    ub(:,2) = uout;
end  

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    nspecies = 1;
    fb = sym(zeros(nspecies+2, 2));

    uin = uinflow(u, mu, eta);
    fin = flux(uin, q, w, v, x, t, mu, eta);
    fb(:,1) = fin *n(1) + tau(1)*(uin-uhat); 

    uout = uoutflow(u, mu, eta);
    fout = flux(uout, q, w, v, x, t, mu, eta);
    fb(:,2) = fout*n(1) + tau(1)*(uout-uhat);
end

function u0 = initu(x, mu, eta)
    nspecies = 1;
    u0 = uinflow(mu,mu,eta);
%     u0 = uoutflow(mu,mu,eta);
end

function uin = uinflow(u, mu, eta)
    nspecies = 1;

    rho = mu(1);
    rhou = mu(2);
    rhoE = mu(3);
    uin = [rho; rhou; rhoE];
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
