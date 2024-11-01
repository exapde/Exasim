function [f,f_udg] = flux(p,udg,param,time)
    %FLUX Volume flux function
    %   [f,fu,fq] = flux(p,u,q,param)
    %
    %      P(N,ND)              Coordinates for N points
    %      U(N,NC)              Unknown vector for N points with NC components
    %      Q(N,NC,ND)           Flux vector for N points with NC components in the
    %                           coordinate directions
    %      PARAM                Parameter list
    %      F(N,NC,ND):          Volume flux at N points
    %      FU(N,NC,ND,NC):      Jacobian of the flux flux vector w.r.t. U
    %      FQ(N,NC,ND,NC,ND):   Jacobian of the flux flux vector w.r.t. Q
    
    [ng,nc] = size(udg);
    nch = 8;
    ns = 5;
    
    gam  = param{1};
    gam1 = gam - 1.0;
                                 
    
    % Nondimensional params
    rho_scale   = param{1};
    u_scale     = param{2};
    rhoe_scale  = param{3};
    T_scale     = param{4};
    mu_scale    = param{5};
    kappa_scale = param{6};
    cp_scale    = param{7};
    L_scale     = param{8};
    Ec = param{9};
    
    
    rmin = 0.0;
    % alphaClip = mu(19);
    % % Conservative Variables
    % for i = 1:ns
    %     rho_i(i) = rmin + lmax(u(i)-rmin,alphaClip); %subspecies density
    %     rho = rho + rho_i(i); %total mixture density
    % end
    rho_i = udg(:,1:ns);
    rho   = sum(rho_i,2);
    rhou  = udg(:,ns+1);
    rhov  = udg(:,ns+2);
    rhoE  = udg(:,ns+3);
    
    drho_dx_i = udg(:,nch + (1:ns));
    drhou_dx  = udg(:,nch + ns+1);
    drhov_dx  = udg(:,nch + ns+2);
    drhoE_dx  = udg(:,nch + ns+2+1);
    drho_dy_i = udg(:,nch + (nch+1:nch+ns));
    drhou_dy  = udg(:,nch + nch+ns+1);
    drhov_dy  = udg(:,nch + nch+ns+2);
    drhoE_dy  = udg(:,nch + nch+ns+2+1);
    av = p(:,3);
    
    rhoinv = 1.0 ./ rho;
    uv = rhou .* rhoinv; %velocity
    vv = rhov .* rhoinv;
    E = rhoE .* rhoinv; %energy
    uTu2   = 0.5*(uv.*uv+vv.*vv);
    
    % TODO: PACK INPUTS
    rho_i_dim = rho_i * rho_scale;
    rho_dim = sum(rho_i_dim,2);
    rhou_dim = rhou * (rho_scale * u_scale);
    rhov_dim = rhov * (rho_scale * u_scale);
    rhoE_dim = rhoE * rhoe_scale;
    rhoe_dim = (rhoE_dim - 0.5 * (rhou_dim.*rhou_dim + rhov_dim.*rhov_dim) ./ rho_dim);
    % TODO: PASS TO MUTATION
    Uin = [rho_i_dim, rhoe_dim];
    [p, T, dT_dr, dT_dre] = mppFlux_rre(ng, 2, Uin);
    
    denom = 1 ./ dT_dre;   
    dT_drho_i = ((uTu2*u_scale^2 ./denom) + dT_dr);
    dT_drhou = (-uv ./ denom * u_scale);
    dT_drhov = (-vv ./ denom * u_scale);
    dT_drhoE = dT_dre;
    
    % TODO: NONDIM pressure_mat
    p = p/rhoe_scale;
    % TODO: GET derivatives of p.
    H = E + p.*rhoinv; %enthalpy
    
    
    % Fluxes
    f = zeros(ng,nch,2);
    for i = 1:ns
        f(:,i,1) = rho_i(:,i) .* uv     + av.*drho_dx_i(:,i);
    end
    f(:,ns + 1,1) = rhou .* uv + p    + av.*drhou_dx;
    f(:,ns + 2,1) = rhov .* uv        + av.*drhov_dx;
    f(:,ns + 3,1) = rhou .* H  + av.*drhoE_dx;
    
    for i = 1:ns
        f(:,i,2) = rho_i(:,i) .* vv     + av.*drho_dy_i(:,i);
    end
    f(:,ns + 1,2) = rhou .* vv        + av.*drhou_dy;
    f(:,ns + 2,2) = rhov .* vv + p    + av.*drhov_dy;
    f(:,ns + 3,2) = rhov .* H         + av.*drhoE_dy;
    
    % f(:,:,1) = [ru+av.*rx, ru.*uv+p+av.*rux, rv.*uv+av.*rvx,   ru.*h+av.*(rEx)];
    % f(:,:,2) = [rv+av.*ry, ru.*vv+av.*ruy,   rv.*vv+p+av.*rvy, rv.*h+av.*(rEy)];
    % f(:,:,1) = [ru, ru.*uv+p, rv.*uv,   ru.*h];
    % f(:,:,2) = [rv, ru.*vv,   rv.*vv+p, rv.*h];
    
    
    if nargout>1
        f_u = zeros(ng,nch,2,nch);
        dp_drho_i = zeros(ng, ns);
        Mw = [14.0067...
        15.9994...
        30.0061...
        28.0134...
        31.9988]; % mix.speciesMw
        RU = 8.314471468617452;
    
        Mw = Mw / 1000.0;
        % NOTE: MAKE SURE THIS IS CALLED ON DIMENSIONAL TEMPERATURE OUTPUTS
        % for i = 1:ns
        %     dp_drho_i(:,i) = (T .* RU./Mw(i)  + pressure_mat(dT_drho_i(:,i), rho_i_dim, Mw)) / rhoe_scale * rho_scale;
        % end
        dp_drho_i = (T .* RU./Mw  + pressure_mat(dT_drho_i, rho_i_dim, Mw)) / rhoe_scale * rho_scale;
        dp_drhou = pressure_mat(dT_drhou, rho_i_dim, Mw) / rhoe_scale * rho_scale * u_scale;
        dp_drhov = pressure_mat(dT_drhov, rho_i_dim, Mw) / rhoe_scale * rho_scale * u_scale;
        dp_drhoE = pressure_mat(dT_drhoE, rho_i_dim, Mw) / rhoe_scale * rhoe_scale;
    
        % derivatives of F_x(1:ns) rho_i u
        for i = 1:ns
            for j = 1:ns
                if j == i
                    f_u(:,i,1,j) = uv .* (rho - rho_i(:,i)) ./ rho;
                else
                    f_u(:,i,1,j) = - rho_i(:,i) .* uv ./ rho;
                end
            end
            % end
            f_u(:,i,1,ns+1) = rho_i(:,i)./rho;
            f_u(:,i,1,ns+2) = zeros(ng,1);
            f_u(:,i,1,ns+3) = zeros(ng,1);
        end
    
        % derivatives of F_x(ns+1) ru^2 + p
        for i = 1:ns
            f_u(:,ns+1,1,i) = -uv.^2 + dp_drho_i(:,i);
        end
        f_u(:,ns+1,1,ns+1) = 2*uv + dp_drhou;
        f_u(:,ns+1,1,ns+2) = dp_drhov;
        f_u(:,ns+1,1,ns+3) = dp_drhoE;
    
        % derivatives of F_x(ns+2) ruv
        for i = 1:ns
            f_u(:,ns+2,1,i) = -uv .* vv;
        end
        f_u(:,ns+2,1,ns+1) = vv;
        f_u(:,ns+2,1,ns+2) = uv;
        f_u(:,ns+2,1,ns+3) = zeros(ng,1);
    
        % derivatives of F_x(ns+3) rhou * H
        for i = 1:ns
            f_u(:,ns+3,1,i) = uv .* dp_drho_i(:,i) - uv .* H;
        end
        f_u(:,ns+3,1,ns+1) = uv .* dp_drhou + H;
        f_u(:,ns+3,1,ns+2) = uv .* dp_drhov;
        f_u(:,ns+3,1,ns+3) = uv .* (1 + dp_drhoE);
    
    
        % derivatives of F_y(1:ns) rho_i v
        for i = 1:ns
            for j = 1:ns
                if j == i
                    f_u(:,i,2,j) = vv .* (rho - rho_i(:,i)) ./ rho;
                else
                    f_u(:,i,2,j) = - rho_i(:,i) .* vv ./ rho;
                end
            end
            % end
            f_u(:,i,2,ns+1) = zeros(ng,1);
            f_u(:,i,2,ns+2) = rho_i(:,i)./rho;
            f_u(:,i,2,ns+3) = zeros(ng,1);
        end
    
        % derivatives of F_y(ns+1) rho uv
        for i = 1:ns
            f_u(:,ns+1,2,i) = -uv .* vv;
        end
        f_u(:,ns+1,2,ns+1) = vv;
        f_u(:,ns+1,2,ns+2) = uv;
        f_u(:,ns+1,2,ns+3) = zeros(ng,1);
    
    
        % derivatives of F_y(ns+2) rv^2 + p
        for i = 1:ns
            f_u(:,ns+2,2,i) = -vv.^2 + dp_drho_i(:,i);
        end
        f_u(:,ns+2,2,ns+1) = dp_drhou;
        f_u(:,ns+2,2,ns+2) = 2*vv + dp_drhov;
        f_u(:,ns+2,2,ns+3) = dp_drhoE;
    
        % derivatives of F_y(ns+3) rhov * H
        for i = 1:ns
            f_u(:,ns+3,2,i) = vv .* dp_drho_i(:,i) - vv .* H;
        end
        f_u(:,ns+3,2,ns+1) = vv .* dp_drhou;
        f_u(:,ns+3,2,ns+2) = vv .* dp_drhov + H;
        f_u(:,ns+3,2,ns+3) = vv .* (1 + dp_drhoE);
    
        % f_u(:,:,1,1) = -[zeros(ng,1), 0.5*((3-gam)*uv.*uv-gam1*vv.*vv), uv.*vv, gam*E.*uv-2*gam1*uv.*af];
        % f_u(:,:,1,2) = -[-ones(ng,1), (gam-3)*uv, -vv, -gam*E+0.5*gam1*(3*uv.*uv+vv.*vv)];
        % f_u(:,:,1,3) = -[zeros(ng,1), gam1*vv, -uv, gam1*uv.*vv];
        % f_u(:,:,1,4) = -[zeros(ng,1), -gam1*ones(ng,1), zeros(ng,1), -gam*uv];    
        
        % f_u(:,:,2,1) = -[zeros(ng,1), uv.*vv, 0.5*((3-gam)*vv.*vv-gam1*uv.*uv), gam*E.*vv-2*gam1*vv.*af];
        % f_u(:,:,2,2) = -[zeros(ng,1), -vv, gam1*uv, gam1*uv.*vv];
        % f_u(:,:,2,3) = -[-ones(ng,1), -uv, (gam-3)*vv,  -gam*E+0.5*gam1*(3*vv.*vv+uv.*uv) ];
        % f_u(:,:,2,4) = -[zeros(ng,1), zeros(ng,1), -gam1*ones(ng,1), -gam*vv];
    
        %f_udg = f_u;
        
        f_q = zeros(ng,nch,2,2*nch);
        for i = 1:nch
            f_q(:,i,1,i) = av;
            f_q(:,i,2,nch+i) = av;
        end    
            
        f_udg = cat(4,f_u,f_q);
    end
    
    