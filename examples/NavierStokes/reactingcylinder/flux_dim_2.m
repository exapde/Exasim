function [f,f_udg] = flux(p,udg,param,time,wdg_1)
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
        
    f_udg = 0;
    [ng,nc] = size(udg);
    nch = 8;
    ns = 5;
    
    gam  = param{1};
    gam1 = gam - 1.0;
                                 
    
    % Nondimensional params
    for i = 1:12; param{i} = 1.0; end
    rho_scale   = 1;
    u_scale     = 1;
    rhoe_scale  = 1;
    T_scale     = 1;
    mu_scale    = 1;
    kappa_scale = 1;
    cp_scale    = 1;
    L_scale     = 1;
    
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
    
    drho_dx_i = -udg(:,nch + (1:ns));
    drhou_dx  = -udg(:,nch + ns+1);
    drhov_dx  = -udg(:,nch + ns+2);
    drhoE_dx  = -udg(:,nch + ns+2+1);
    drho_dy_i = -udg(:,nch + (nch+1:nch+ns));
    drhou_dy  = -udg(:,nch + nch+ns+1);
    drhov_dy  = -udg(:,nch + nch+ns+2);
    drhoE_dy  = -udg(:,nch + nch+ns+2+1);
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
    rhoe_dim = rhoE_dim - 0.5 * (rhou_dim.*rhou_dim + rhov_dim.*rhov_dim) ./ rho_dim;
    % TODO: PASS TO MUTATION
    Uin = [rho_i_dim, rhoe_dim];
    [wdg, dwdg] = state2mutation_my_temp(udg, param);
    if nargin <5
        T = wdg(:,1);
    else
        T = wdg_1(:,1);
    end
    T_dim = T*T_scale;
    p = wdg(:,2);
    p_dim = p*rhoe_scale;
    
    dT_drho_i = dwdg(:,1,1:5);
    dT_drhou  = dwdg(:,1,6);
    dT_drhov  = dwdg(:,1,7);
    dT_drhoE  = dwdg(:,1,8);
    
    % TODO: NONDIM PRESSURE
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
    
    
    % Viscous fluxes
    beta = 0;
    Ec          = 1.0;
    Pr          = 1.0;
    Re          = 1.0;
    [species_thermo_structs, Mw, RU] = thermodynamicsModels();
    [blottner_structs, gupta_structs, gupta_mu_structs, gupta_kappa_structs] = transport();

    drho_dx = sum(drho_dx_i,2);
    drho_dy = sum(drho_dy_i,2);
    du_dx = (drhou_dx - drho_dx.*uv).*rhoinv;
    dv_dx = (drhov_dx - drho_dx.*vv).*rhoinv;
    du_dy = (drhou_dy - drho_dy.*uv).*rhoinv;
    dv_dy = (drhov_dy - drho_dy.*vv).*rhoinv;
    uTu2      = 0.5 * (uv .* uv + vv .* vv);
    duTu2_dx  = uv .* du_dx + vv .* dv_dx; 
    duTu2_dy  = uv .* du_dy + vv .* dv_dy;
    E = rhoE .* rhoinv; 
    uv = rhou .* rhoinv; %velocity
    vv = rhov .* rhoinv;

    dre_drho  = -uTu2;
    dre_duTu2 = -rho;
    dre_drhoE = 1.0*ones(size(p,1));

    for ig = 1:size(p,1)

        X = X_i(rho_i_dim(ig,:)',Mw);

        Y = Y_i(rho_i_dim(ig,:)');
        denom = sum(rho_i_dim(ig,:)') * mixtureFrozenCvMass(T_dim(ig), Mw, Y, species_thermo_structs);
        e_i = getEnergiesMass(T_dim(ig), Mw, species_thermo_structs);
        dT_drho_i_dim = -e_i(:) ./ denom;
        dT_drhoe_dim = 1.0 / denom;

        dT_drho_i = dT_drho_i_dim / T_scale * rho_scale;
        dT_drhoe = dT_drhoe_dim / T_scale * rhoe_scale;

        dre_dx    = dre_drho(ig) * drho_dx(ig) + dre_duTu2(ig) * duTu2_dx(ig) + dre_drhoE(ig) * drhoE_dx(ig);
        dre_dy    = dre_drho(ig) * drho_dy(ig) + dre_duTu2(ig) * duTu2_dy(ig) + dre_drhoE(ig) * drhoE_dy(ig);
        dT_dx     = sum(dT_drho_i .* drho_dx_i(ig,:)') +  dT_drhoe * dre_dx;
        dT_dy     = sum(dT_drho_i .* drho_dy_i(ig,:)') +  dT_drhoe * dre_dy;

        D_vec = averageDiffusionCoeffs(T_dim(ig), X, Y, Mw, p_dim(ig), gupta_structs);
        h_vec = getEnthalpiesMass(T_dim(ig), Mw, species_thermo_structs);

        mu_i = speciesViscosities(T_dim(ig), gupta_mu_structs);
        % lambda_i = speciesConductivities(T, gupta_kappa_structs);
        phi_i = euckenPhi(mu_i, Mw, X);
        mu_d_dim = wilkeMixture(mu_i, X, phi_i);
        lambda_i = mu_d_dim * 3/2 .* getCpsMass(T_dim(ig), Mw, species_thermo_structs);
        kappa_dim = wilkeMixture(lambda_i, X, phi_i);

        h_scale = u_scale^2;
        D_scale = u_scale;

        mu_d = mu_d_dim / mu_scale;
        kappa = kappa_dim / kappa_scale;
        D_vec = D_vec / D_scale;
        h_vec = h_vec / h_scale;

        %%%%%%%% Calculation of J_i
        dY_dx_i = (drho_dx_i(ig,:)' * rho(ig) - rho_i(ig,:)' * drho_dx(ig)) * rhoinv(ig) * rhoinv(ig);
        dY_dy_i = (drho_dy_i(ig,:)' * rho(ig) - rho_i(ig,:)' * drho_dy(ig)) * rhoinv(ig) * rhoinv(ig);

        J_i_x = -rho(ig) * D_vec .* dY_dx_i + rho_i(ig,:)' .* sum(D_vec .* dY_dx_i);
        J_i_y = -rho(ig) * D_vec .* dY_dy_i + rho_i(ig,:)' .* sum(D_vec .* dY_dy_i);

        %%%%%%%% Stress tensor tau
        txx = mu_d * 2.0/3.0 * (2 * du_dx(ig) - dv_dy(ig)) / Re ;
        txy = mu_d * (du_dy(ig) + dv_dx(ig)) / Re;
        tyy = mu_d * 2.0/3.0 * (2 * dv_dy(ig) - du_dx(ig)) / Re;

        % VISCOUS FLUX
        
        for i = 1:ns
            fv(ig,i,1) = -J_i_x(i); 
            fv(ig,i,2) = -J_i_y(i);
        end

        fv(ig,ns + 1, 1) = txx;
        fv(ig,ns + 2, 1) = txy;
        fv(ig,ns + 3,1) = (uv(ig) * txx + vv(ig) * txy) - (sum(h_vec.*J_i_x) - kappa *dT_dx / (Re*Pr*Ec));
        
        fv(ig,ns + 1, 2) = txy;
        fv(ig,ns + 2, 2) = tyy;
        fv(ig,ns + 3,2) = (uv(ig) * txy + vv(ig) * tyy) - (sum(h_vec.*J_i_y) - kappa *dT_dy/ (Re*Pr*Ec));
    end
    
    f = f- fv;
end