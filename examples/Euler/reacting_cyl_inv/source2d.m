function [sr,sr_udg] = source(p,udg,param,time)

    [ng,nc] = size(udg);
    nch = 8;
    ns = 5;
    
    sr = zeros(ng,nch);
    sr_udg = zeros(ng,nch,nc); 
    
    rho_scale   = param{1};
    u_scale     = param{2};
    rhoe_scale  = param{3};
    T_scale     = param{4};
    mu_scale    = param{5};
    kappa_scale = param{6};
    cp_scale    = param{7};
    L_scale     = param{8};
    omega_scale = rho_scale * u_scale / L_scale;
    
                                                 
    r_i  = udg(:,1:ns);
    ru   = udg(:,ns+1);
    rv   = udg(:,ns+2);
    rE   = udg(:,ns+3);
    r = sum(r_i, 2);
    
    r1   = 1./r;
    uv   = ru.*r1;
    vv   = rv.*r1;
    E    = rE.*r1;
    uTu2   = 0.5*(uv.*uv+vv.*vv);
    
    % TODO: map to rho_i, rhoe
    % TODO: dimensionalize 
    rho_i_dim = r_i * rho_scale;
    rho_dim = sum(rho_i_dim, 2);
    rhou_dim = ru * (rho_scale * u_scale);
    rhov_dim = rv * (rho_scale * u_scale);
    rhoE_dim = rE * rhoe_scale;
    rhoe_dim = rhoE_dim - 0.5 * (rhou_dim.*rhou_dim + rhov_dim.*rhov_dim) ./ rho_dim;
    
    Uin = [rho_i_dim, rhoe_dim];
    [w_i, dw_dr, dw_dT, dT_dr, dT_dre] = mppSource(ng,1e-6,Uin);
    
    % TODO: THIS COULD BE WRONG BUT IT'S A BIT SIMPLER...NEED TO BE CAREFUL HERE 
    denom = 1 ./ dT_dre;
    dT_drho_i = (uTu2.*u_scale^2./denom) + dT_dr;
    dT_drhou = -uv ./ denom * u_scale;
    dT_drhov = -vv ./ denom * u_scale;
    dT_drhoE = dT_dre;
    
    sr(:,1:ns) = w_i / omega_scale;
    
    if nargout>1
        f_u = zeros(ng,nch,nch);
        % Note: dimensional terms here
        for i = 1:ns
            for j = 1:ns
                dw_drho(:,i,j) = dw_dr(:,i,j) + dw_dT(:,i) .* dT_drho_i(:,j);
            end
        end
        % dw_drho = dw_dr + dw_dT .* dT_drho_i;
        dw_drhou = dw_dT .* dT_drhou;
        dw_drhov = dw_dT .* dT_drhov;
        dw_drhoE = dw_dT .* dT_drhoE;
        %
    
        % f_u(:,1:ns,1:ns) = permute(dw_drho, [1 3 2]) * rho_scale / omega_scale;
        f_u(:,1:ns,1:ns) = dw_drho * rho_scale / omega_scale;
    
        % f_u(:,1:ns,2) = dw_drho(:,2) * rho_scale / omega_scale;
        % f_u(:,1:ns,3) = dw_drho(:,3) * rho_scale / omega_scale;
        % f_u(:,1:ns,4) = dw_drho(:,4) * rho_scale / omega_scale;
        % f_u(:,1:ns,5) = dw_drho(:,5) * rho_scale / omega_scale;
    
        f_u(:,1:ns,ns+1) = dw_drhou * (rho_scale * u_scale) / omega_scale;
        f_u(:,1:ns,ns+2) = dw_drhov * (rho_scale * u_scale) / omega_scale;
        f_u(:,1:ns,ns+3) = dw_drhoE * rhoe_scale / omega_scale;
        
        f_q = zeros(ng,nch,2*nch);
    
        
        sr_udg = cat(3,f_u,f_q);
    end
    