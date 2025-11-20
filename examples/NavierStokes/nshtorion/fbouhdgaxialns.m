function fb = fbouhdgaxialns(u, q, w, v, x, t, mu, eta, uhat, n, tau)

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

    % iso-thermal wall boundary condition    
    f_iso = 0*u;
    f_iso(1) = u(1) - uhat(1); % extrapolate density
    f_iso(2) = 0.0  - uhat(2); % zero velocity
    f_iso(3) = 0.0  - uhat(3); % zero velocity           
    f_iso(4) = -uhat(4) + uhat(1)*TisoW; % set temperature to Twall
      
    % iso-thermal wall boundary condition    
    f3 = 0*u;
    f3(1) = u(1) - uhat(1); % extrapolate density
    f3(2) = 0.0  - uhat(2); % zero velocity
    f3(3) = 0.0  - uhat(3); % zero velocity           
    f3(4) = -uhat(4) + uhat(1)*v(2);

    % symmetry boundary condition    
    r = uhat(1);
    ru = uhat(2);
    rv = uhat(3);
    rE = uhat(4);
    r1 = 1/r;
    uv = ru*r1;        
    vv = rv*r1;    
    ke = 0.5*(uv*uv+vv*vv);
    p = gam1*(rE-r*ke);

    ry = q(5);
    ruy = q(6);    
    rvy = q(7);
    rEy = q(8);    
    uy = (ruy - ry*uv)*r1;    
    vy = (rvy - ry*vv)*r1;
    qy = uv*uy + vv*vy;
    py = gam1*(rEy - ry*ke - r*qy);
    Ty = 1/gam1*(py*r - p*ry)*r1^2;

    f_sym = 0*u;
    f_sym(1) = ry + tau*(u(1) - uhat(1));  % dr/dy = 0
    f_sym(2) = uy + tau*(u(2) - uhat(2));  % du/dy = 0
    f_sym(3) = vy + tau*(0.0 - uhat(3));   % v = 0
    f_sym(4) = Ty + tau*(u(4) - uhat(4));  % dT/dy = 0
    
    % slip wall condition                
    ru = u(2);
    rv = u(3);
    nx = n(1);
    ny = n(2);   
    run = ru*nx + rv*ny;        
    uinf = u;
    uinf(2) = uinf(2) - nx.*run;
    uinf(3) = uinf(3) - ny.*run;
    f_slip = tau*(uinf - uhat);    
    
    % zero gradient condition
    q = q(:);
    f_grad = q(1:4)*nx + q(5:8)*ny + tau*(u(:) - uhat(:));
        
    % supersonic inflow, supersonic outflow, isothermal, symmetry, slip wall, gradient                 
    fb = [f_in f_out f_iso f_grad];
end

