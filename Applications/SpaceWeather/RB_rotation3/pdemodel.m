function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym([1.0; 1.0; 1.0; 1.0]); 
end

function f = flux(u, q, w, v, x, t, mu, eta)
    gam = mu(1);
    gam1 = gam - 1.0;
    Re = mu(2);
    Pr = mu(3);
    Minf = mu(4);
    
    cp = mu(7);
    cv = cp/gam;
    R = cp-cv;
    
%     muRef = 1/Re;
    muRef = mu(15);
    c23 = 2.0/3.0;
    
    Tref = 273;
    
    r = u(1);
    ru = u(2);
    rv = u(3);
    rE = u(4);
    
    rx = q(1);
    rux = q(2);
    rvx = q(3);
    rEx = q(4);
    ry = q(5);
    ruy = q(6);
    rvy = q(7);
    rEy = q(8);
    
    r1 = 1.0/r;
    uv = ru*r1;
    vv = rv*r1;
    E = rE*r1;
    ke = 0.5*(uv*uv+vv*vv);
    p = gam1*(rE-r*ke);
    h = E+p*r1;
    
    % Inviscid fluxes
    fi = [ru, ru*uv+p, rv*uv, ru*h, ...
            rv, ru*vv, rv*vv+p, rv*h];
        
    ux = (rux - rx*uv)*r1;
    vx = (rvx - rx*vv)*r1;
    kex = uv*ux + vv*vx;
    px = gam1*(rEx - rx*ke - r*kex);
    Tx = (px*r - p*rx)*r1^2/R;
    uy = (ruy - ry*uv)*r1;
    vy = (rvy - ry*vv)*r1;
    key = uv*uy + vv*vy;
    py = gam1*(rEy - ry*ke - r*key);
    Ty = (py*r - p*ry)*r1^2/R;

    % Viscosity
    T = p/(R*r);
%     muSuth = getViscosity(muRef,Tref,T,1);
    muSuth = muRef*sqrt(T/R);
    fc = muSuth*cp/Pr;

    %Viscous fluxes
    txx = muSuth*c23*(2*ux - vy);
    txy = muSuth*(uy + vx);
    tyy = muSuth*c23*(2*vy - ux);
    fv = [0, txx, txy, uv*txx + vv*txy + fc*Tx, ...
          0, txy, tyy, uv*txy + vv*tyy + fc*Ty];
    f = fi+fv;
    f = reshape(f,[4,2]);    
end

function s = source(u, q, w, v, x, t, mu, eta)
    factor = eta(end-1);
    t0 = eta(end);
    f_t = 0.5*tanh(factor*(t-t0))+0.5;
    s_EUV = f_t*EUVsource(u, x, t, mu, eta);
    
    r = sqrt(x(1)^2 + x(2)^2);
    sca = (x(1)*u(2) + x(2)*u(3))/r;
    
    R0 = mu(12);
    gravity0 = mu(5);
    gravity = gravity0*(R0/r)^2;
    
    omega = mu(6);
    s = [sym(0.0); -gravity*u(1)*x(1)/r + 2.0*omega*u(3) + x(1)*omega^2; ...
                   -gravity*u(1)*x(2)/r - 2.0*omega*u(2) + x(2)*omega^2; ...
                   -gravity*sca + omega^2*sca*r + s_EUV];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    fi = f(:,1)*n(1) + f(:,2)*n(2) + tau*(u-uhat); % numerical flux at freestream boundary
    
    fa = fi;
    fa(4) = 0.0;
    
    fw = fi;
    fw(1) = 0.0;   % zero velocity 
    fb = [fw fa];  % isothermal, "zero gradient"
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)

    gam = mu(1);
    gam1 = gam - 1.0;
    cp = mu(7);
    cv = cp/gam;
    R = cp-cv;
    gravity0 = mu(5);
    pbot = mu(9);
    Tbot = mu(10);
    Ttop = mu(11);
    
    Rbot = mu(12);
    Rtop = mu(13);
    
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
    Temp = p/(R*r);
    
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
    gravity = gravity0*Rbot^2/Rtop^2;
    ptop = pbot*exp(-20*gravity/(R*Ttop)*log(1 + Ttop/Tbot*(exp(0.05*(Rtop-Rbot)-1))));
    rtop = ptop/(R*Ttop);
    
    uinf = sym([rtop; 0.0; 0.0; rtop*cv*Temp]); % freestream flow
    uinf = uinf(:);
    u = u(:);          % state variables 
    ui = 0.5*((u+uinf) + An*(u-uinf));  % Riemann solution

    % Isothermal Wall?
    utw1 = u;
    utw1(2:3) = 0.0;
    utw1(4) = u(1)*cv*Tbot;

    ub = [utw1 ui]; % isothermal, "zero gradient"
end


function u0 = initu(x, mu, eta)
    gam = mu(1);
    gravity0 = mu(5);
    cp = mu(7);
    cv = cp/gam;
    R = cp-cv;

    pbot = mu(9);
    Tbot = mu(10);
    Ttop = mu(11);
    R0 = mu(12);
    
    r = sqrt(x(1)^2 + x(2)^2);
    
    gravity = gravity0*R0^2/r^2;

    A = Ttop;
    B = -(Ttop-Tbot)*exp(0.05*R0);
    C = 0.05;
    f = log(A*exp(C*r)+B)/(A*C);
    f0 = log(A*exp(C*R0)+B)/(A*C);

    T = A+B*exp(-C*r);
    p = pbot*exp(-gravity*(f-f0)/R);
    rho = p/(R*T);

    u0 = [ rho, 0.0, 0.0, p/(gam-1.0)];    
end
