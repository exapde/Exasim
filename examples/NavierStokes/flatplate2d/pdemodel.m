function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
pde.fbouhdg = @fbouhdg;
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
    Re1 = 1/Re;
    M2 = Minf^2;
    c23 = 2.0/3.0;
    fc = 1/(gam1*M2*Re*Pr);
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
    r1 = 1/r;
    uv = ru*r1;
    vv = rv*r1;
    E = rE*r1;
    ke = 0.5*(uv*uv+vv*vv);
    p = gam1*(rE-r*ke);
    h = E+p*r1;
    fi = [ru, ru*uv+p, rv*uv, ru*h, ...
            rv, ru*vv, rv*vv+p, rv*h];
    ux = (rux - rx*uv)*r1;
    vx = (rvx - rx*vv)*r1;
    kex = uv*ux + vv*vx;
    px = gam1*(rEx - rx*ke - r*kex);
    Tx = gam*M2*(px*r - p*rx)*r1^2;
    uy = (ruy - ry*uv)*r1;
    vy = (rvy - ry*vv)*r1;
    key = uv*uy + vv*vy;
    py = gam1*(rEy - ry*ke - r*key);
    Ty = gam*M2*(py*r - p*ry)*r1^2;
    txx = Re1*c23*(2*ux - vy);
    txy = Re1*(uy + vx);
    tyy = Re1*c23*(2*vy - ux);
    fv = [0, txx, txy, uv*txx + vv*txy + fc*Tx, ...
          0, txy, tyy, uv*txy + vv*tyy + fc*Ty];
    f = fi+fv;
    f = reshape(f,[4,2]);      
end

function s = source(u, q, w, v, x, t, mu, eta)

s = [sym(0.0); sym(0.0); sym(0.0); sym(0.0)];    
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    fb = 0*uhat; 
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = 0*uhat; 
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)

    gam = mu(1);
    gam1 = gam - 1.0;
    r = uhat(1);
    ru = uhat(2);
    rv = uhat(3);
    rE = uhat(4);
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
    Lambda = [ tanh(1e2*(un-a)) , 0 , 0 , 0 ;...
                     0 , tanh(1e2*(un)) , 0 , 0 ;...
                     0 , 0 , tanh(1e2*(un)) , 0 ;...
                     0 , 0 , 0 , tanh(1e2*(un+a)) ];
%     E = (K * Lambda * Kinv);
%     An = (Tinv * E * T);
    L = simplify(Tinv * K);
    R = simplify(Kinv * T);
    An = simplify(L * Lambda * R);

    % freestream boundary condition
    uinf = sym(mu(5:8)); % freestream flow variables
    uinf = uinf(:);
    u = u(:);          % state variables 
    fb1 = 0.5*((u+uinf) + An*(u-uinf)) - uhat;     
    
    % wall boundary condition    
    fb2 = 0*fb1;
    fb2(1) = u(1) - uhat(1); % extrapolate density
    fb2(2) = 0.0  - uhat(2); % zero velocity
    fb2(3) = 0.0  - uhat(3); % zero velocity           
    f = flux(uhat, q, w, v, x, t, mu, eta);
    fb2(4) = f(4,1)*n(1) + f(4,2)*n(2) + tau*(u(4)-uhat(4)); % zero heat flux
       
    % outflow boundary condition  
    Minf = mu(4);
    pinf = 1/(gam*Minf^2); 
    uinf = [u(1), u(2), u(3), pinf/(gam-1) + 0.5*u(2)*u(2)/u(1) + 0.5*u(3)*u(3)/u(1)];
    uinf = uinf(:);
    fb3 = 0.5*((u+uinf) + An*(u-uinf)) - uhat;     
    
    fb = [fb1 fb2 fb3];    
end

function u0 = initu(x, mu, eta)
    u0 = sym(mu(5:8)); % freestream flow   
end




