function fb = fbouhdgaxialns(u, q, w, v, x, t, mu, eta, uhat, n, tau)

    % symmetry boundary condition    
    r = uhat(1);
    ru = uhat(2);
    rv = uhat(3);
    rE = uhat(4);
    ry = q(5);
    ruy = q(6);    
    rvy = q(7);
    rEy = q(8);  

    % regularization mueters
    alpha = 4.0e3;
    rmin = 5.0e-2;
    pmin = 2.0e-3;

    % Regularization of rho (cannot be smaller than rmin)
    r = rmin + lmax(r-rmin,alpha);
    % Density sensor
    dr = atan(alpha*(r - rmin))/pi + (alpha*(r - rmin))/(pi*(alpha^2*(r - rmin)^2 + 1)) + 1/2;
    r1 = 1/r;
    uv = ru*r1;
    vv = rv*r1;
    E = rE*r1;

    gam = 1.4;        % Thermodynamic quantities are determined consistently from the equation of state and used to couple the energy to the momentum equation in the following manner: energy->temperature->pressure
    q = 0.5*(uv*uv+vv*vv);
    e = E-q;
    T = e*(gam-1);
    p = r*T;
    H = E + T;

    p = pmin + lmax(p-pmin,alpha);

    uy = (ruy - ry*uv)*r1;      % Product rule -- we need to isolate the derivative of velocity from the derivative of momentum
    vy = (rvy - ry*vv)*r1;
    Ey = (rEy - ry*E)*r1;
    Ty = (gam-1)*(Ey - uv*uy - vv*vy);

    f_sym = 0*u;
    % Two choices for symmetry BC: 0 velocity or 0 velocity gradient.
    % This block is for zero velocity gradient
    f_sym(1) = ry + tau*(u(1) - uhat(1));  % dr/dy = 0
    f_sym(2) = uy + tau*(u(2) - uhat(2));  % du/dy = 0
    f_sym(3) = vy + tau*(u(3) - uhat(3));  % dv/dy = 0
    f_sym(4) = Ty + tau*(u(4) - uhat(4));  % dT/dy = 0

    % % Zero velocity at the symmetry condition: other possible boundary conditions to try
    % f_sym(1) = ry + tau*(u(1) - uhat(1));  % dr/dy = 0
    % f_sym(2) = 0.0  - nx*uhat(2); % zero velocity component in x
    % f_sym(3) = 0.0  - ny*uhat(3); % zero velocity component in y
    % f_sym(4) = Ty + tau*(u(4) - uhat(4));  % dT/dy = 0

    % Other possible boundary conditions to try
    % f_sym(1) = ry + tau*(u(1) - uhat(1));  % dr/dy = 0
    % f_sym(2) = uy + tau*(u(2) - uhat(2));  % du/dy = 0
    % f_sym(3) = vy + tau*(0.0 - uhat(3));   % u = 0
    % f_sym(4) = Ty + tau*(u(4) - uhat(4));  % dT/dy = 0

    % Outflow boundary
    nx = n(1);
    ny = n(2);
    gam1 = gam-1;
    a = sqrt(gam*p*r1);           % Using nondimensional acoustic velocity

    run = ru*nx + rv*ny;
    rut = -ru*ny + rv*nx;
    un = run*r1;
    ut = rut*r1;
    K = [ 1 , 1 , 0 , 1 ;...            % Euler flux jacobian
        un-a , un , 0 , un+a ;...
        ut , ut , 1 , ut ;...
        H - un*a , (1/2)*(un^2 + ut^2) , ut , H+un*a ];
    Kinv = (gam1/(2*a^2))*[ H + (a/gam1)*(un-a) , -(un+a/gam1) , -ut , 1 ;...       % Inverse of Euler flux jacobian
                            -2*H + (4/gam1)*a^2 , 2*un , 2*ut , -2 ;...
                            -2*(ut*a^2)/gam1 , 0 , 2*(a^2)/gam1 , 0 ;...
                            H - a*(un+a)/gam1 , -un+a/gam1 , -ut , 1 ];
    T = [ 1 , 0 , 0 , 0;...     % Rotation matrix for boundary in 2D
        0 , nx , ny , 0;...
        0 , -ny , nx , 0;...
        0 , 0 , 0 , 1];
    Tinv = [ 1 , 0 , 0 , 0;...     % Inverse rotation matrix for boundary in 2D
            0 , nx ,-ny , 0;...
            0 , ny , nx , 0;...
            0 , 0 , 0 , 1];
    Lambda = [ tanh(1e2*(un-a)) , 0 , 0 , 0 ;...        % We don't need to multiply by r because r only scales the eigenvalues and is is always >0
                    0 , tanh(1e2*(un)) , 0 , 0 ;...
                    0 , 0 , tanh(1e2*(un)) , 0 ;...
                    0 , 0 , 0 , tanh(1e2*(un+a)) ];

    L = simplify(Tinv * K);
    R = simplify(Kinv * T);
    An = simplify(L * Lambda * R);

    % Wall boundary condition    
    fwall = 0*f_sym;
    fwall(1) = u(1) - uhat(1); % extrapolate density
    fwall(2) = 0.0  - nx*uhat(2); % zero velocity component in x
    fwall(3) = 0.0  - ny*uhat(3); % zero velocity component in y
    fwall(4) = u(4) - uhat(4); % extrapolate energy      Different from NS

    % Outflow boundary condition  
    pinf=1.225853658536585;     % Outlet pressure = reference pressure by construction of the nondimensionalization
    uinf = [u(1), u(2), u(3), pinf/(gam-1) + 0.5*u(2)*u(2)/u(1) + 0.5*u(3)*u(3)/u(1)];
    uinf = uinf(:);
    foutflow = 0.5*((u(1:4)+uinf(1:4)) + An*(u(1:4)-uinf(1:4))) - uhat(1:4);     

    fb = [f_sym fwall fwall foutflow];
end
