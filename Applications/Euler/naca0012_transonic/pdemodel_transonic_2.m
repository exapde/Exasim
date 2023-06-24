function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
pde.avfield = @avfield;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym([1.0; 1.0; 1.0; 1.0]); 
end

function f = flux(u, q, w, v, x, t, mu, eta)
    f = getfluxav2d_2(u,q,v,mu);
    f = reshape(f,[4,2]);        
end

function f = avfield(u, q, w, v, x, t, mu, eta)
    f = getavfield2d_pb_2(u,q,v,mu);
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
    r = u(1);
    ru = u(2);
    rv = u(3);
    rE = u(4);
    nx = n(1);
    ny = n(2);
    
    r1 = 1/r;
    uv = ru*r1;
    vv = rv*r1;
    E = rE*r1;
    p = gam1*(rE-r*0.5*(uv*uv+vv*vv));
    h = E+p*r1;
    a = sqrt(gam*p*r1);
    
    run = ru*nx + rv*ny;
    rut = -ru*ny + rv*nx;
    un = run/r;
    ut = rut/r;
    
    K = [ 1 , 1 , 0 , 1 ;...
          un-a , un , 0 , un+a ;...
          ut , ut , 1 , ut ;...
          h - un*a , (1/2)*(un^2 + ut^2) , ut , h+un*a ];
    Kinv = (gam1/(2*a^2))*[ h + (a/gam1)*(un-a) , -(un+a/gam1) , -ut , 1 ;...
                            -2*h + (4/gam1)*a^2 , 2*un , 2*ut , -2 ;...
                            -2*(ut*a^2)/gam1 , 0 , 2*(a^2)/gam1 , 0 ;...
                            h - a*(un+a)/gam1 , -un+a/gam1 , -ut , 1 ];
    T = [ 1 , 0 , 0 , 0;...
          0 , nx , ny , 0;...
          0 , -ny , nx , 0;...
          0 , 0 , 0 , 1];
    Tinv = [ 1 , 0 , 0 , 0;...
             0 , nx ,-ny , 0;...
             0 , ny , nx , 0;...
             0 , 0 , 0 , 1];
    Lambda = [ tanh(1e3*(un-a)) , 0 , 0 , 0 ;...
                     0 , tanh(1e3*(un)) , 0 , 0 ;...
                     0 , 0 , tanh(1e3*(un)) , 0 ;...
                     0 , 0 , 0 , tanh(1e3*(un+a)) ];
    E = (K * Lambda * Kinv);
    An = (Tinv * E * T);
    
    % freestream boundary condition
    uinf = sym(mu(3:6)); % freestream flow
    uinf = uinf(:);
    u = u(:);          % state variables 
    
    % freestream
    ui = 0.5*((u+uinf) + An*(u-uinf));  % Riemann solution
    
    % adiabatic wall
%     uaw = u;
%     uaw(2:3) = 0; % zero velocity

    uaw = u;
    uaw(2) = ru - run*nx;
    uaw(3) = rv - run*ny;

    % Isothermal Wall
%     Tinf = mu(9);
%     Tref = mu(10);
%     Twall = mu(11);
%     TisoW = Twall/Tref * Tinf;
    utw = u(:);
    utw(2:3) = 0;
    utw(4) = u(1)*0;
    
    % Slip wall
    usw = u;
    usw(2) = u(2) - nx * (u(2)*nx + u(3)*ny);
    usw(3) = u(3) - ny * (u(2)*nx + u(3)*ny);
    
    % freestream, adiabatic wall, isothermal wall, adiabatic slip wall, supersonic inflow, supersonic outflow
    ub = [ui uaw utw usw uinf u]; 
end

function u0 = initu(x, mu, eta)
    u0 = sym(mu(3:6)); % freestream flow   
end




