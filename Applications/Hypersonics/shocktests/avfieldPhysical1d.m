function avField = avfieldPhysical1d(u, q, w, v, x, t, mu, eta, nspecies, nd)
    % Mutation outputs
    p = eta(1);
    gam = eta(2);
    gam1 = gam - 1;
    avk_b_coeff = eta(3); % placeholder for avk coefficient. Something with units cp / Pr 
    
    % artificial viscosity
    porder = mu(15);

    % regularization parameters for the bulk viscosity
    kb = mu(16);
    sb0   = mu(17); %0.02 
    sbmin = 0.0;
    sbmax = mu(18);% / sqrt(gam*gam - 1.0); %2.5 TODO: unsure that this should be changing 

    % regularization parameters
    alpha = 1.0e12;
    rmin = 0.0;
    Hmin = 1.0e-4;

    % mesh size
    hm = v(4+nspecies);

    % Get base variables
    rvec = u(1:nspecies);
    ru = u(nspecies+1);
    % rv = udg(3);
    rE = u(nspecies+nd+1);

    rx = sum(-q(1:nspecies));
    rux = -q(nspecies+1);

    % Regularization 
    rvec = rmin + lmax(rvec-rmin,alpha); % need to double check regularizatino here
    
    % Simple derived working quantities
    r = sum(rvec);
    r1 = 1./r;
    uv = ru.*r1;
    E = rE * r1; %energy
    H = E + p*r1; %enthalpy
    H = Hmin + lmax(H - Hmin, alpha);

    % Critical speed of Sound
    c_star = sqrt((2.*gam1.*H) ./ (gam+1)); %TODO: this is fine for 1 temp but not 2

    % Computing derivatives for the sensors
    ux = (rux - rx.*uv).*r1;
    div_v = (ux); %TODO: check signs. Actually this probably is important, we want to make sure it's applied in negative v. pos dilitation
                  %      pretty sure my signs are okay. There is something strange about the code I was given. Truly don't understand the signs 
    % limit  divergence and vorticity
    sigm = 1e4;
    div_v = limiting(div_v,-sigm,sigm,alpha,-sigm);

    % Dilatation Sensor sb
    sb = - (hm./porder) .* (div_v./c_star);
    sb = limiting(sb,sbmin,sbmax,alpha,sb0); % TODO: what should sbmin, sbmax, alpha, and sb0 be 

    % Artificial Bulk viscosity
    avb = r * (kb.*hm./(porder)) .* sqrt(uv^2 + c_star^2) .* sb;
    
    % Artificial conductivity
    avk_T = 0.0;                  %  Triggered by temperature sensor
    avk_b =  avk_b_coeff * avb;   %  Triggered by bulk viscosity 
    % Assign artificial viscosities
    avField(1) = avb;             %  bulk
    avField(2) = 0.0;             %  shear
    avField(3) = avk_b + avk_T;   %  thermal
    avField(4:3+nspecies) = 0.0;  %  species
end
