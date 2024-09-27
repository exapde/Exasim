function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
pde.avfield = @avfield;
pde.fbouhdg = @fbouhdg;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym([1.0; 1.0; 1.0; 1.0]); 
end

function f = flux(u, q, w, v, x, t, mu, eta)
    f = getfluxav2d(u,q,v,mu);
    f = reshape(f,[4,2]);        
end

function f = avfield(u, q, w, v, x, t, mu, eta)
    f = getavfield2d(u,q,v,mu);
end

function s = source(u, q, w, v, x, t, mu, eta)
    s = [sym(0.0); sym(0.0); sym(0.0); sym(0.0)];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)

    f = flux(uhat, q, w, v, x, t, mu, eta);
    fi = f(:,1)*n(1) + f(:,2)*n(2) + tau*(u-uhat); % numerical flux at freestream boundary
    
    % adiabatic wall
    faw = fi;
    faw(1) = 0.0;   % zero velocity 
    faw(end) = 0.0; % adiabatic wall -> zero heat flux
    
    % Flux Thermal Wall
    ftw = fi;
    ftw(1) = 0.0;
    
    % freestream, adiabatic wall, isothermal wall, adiabatic slip wall, supersonic inflow, supersonic outflow
    fb = [fi faw ftw faw fi fi]; 
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)

    gam = mu(1);
    gam1 = gam - 1.0;
    
    % freestream boundary condition
    uinf = sym(mu(5:8)); % freestream flow
    uinf = uinf(:);
    u = u(:);          % state variables 

    nx = n(1); ny = n(2);

    % Isothermal Wall
    Tinf = mu(9);
    Tref = mu(10);
    Twall = mu(11);
    TisoW = Twall/Tref * Tinf;
    utw = u(:);
    utw(2:3) = 0;
    utw(4) = u(1)*TisoW;
    
    % Slip wall
    usw = u;
    usw(2) = u(2) - nx * (u(2)*nx + u(3)*ny);
    usw(3) = u(3) - ny * (u(2)*nx + u(3)*ny);
    
    % freestream, adiabatic wall, isothermal wall, adiabatic slip wall, supersonic inflow, supersonic outflow
    ub = [uinf uinf utw usw uinf u]; 
end
function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)


    gam = mu(1);
    gam1 = gam - 1.0;
    Tinf = mu(9);
    Tref = mu(10);
    Twall = mu(11);
    TisoW = Twall/Tref * Tinf;    
    uinf = sym(mu(5:8)); % freestream flow
    uinf = uinf(:);

    f_out = u - uhat;
    f_in = uinf - uhat;

    % wall boundary condition    
    f1 = 0*u;
    f1(1) = u(1) - uhat(1); % extrapolate density
    f1(2) = 0.0  - uhat(2); % zero velocity
    f1(3) = 0.0  - uhat(3); % zero velocity           
    f1(4) = -uhat(4) +uhat(1)*TisoW;
    % f = flux(uhat, q, w, v, x, t, mu, eta);
    % f1(4) = f(4,1)*n(1) + f(4,2)*n(2) + tau*(u(4)-uhat(4)); % zero heat flux
    % freestream, adiabatic wall, isothermal wall, adiabatic slip wall, supersonic inflow, supersonic outflow
    fb = [f_in f_in f1 f1 f_in f_out];
end

function u0 = initu(x, mu, eta)
    u0 = sym(mu(5:8)); % freestream flow   
end




