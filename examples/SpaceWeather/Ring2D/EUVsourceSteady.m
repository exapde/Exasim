function s_EUV = EUVsourceSteady(u, x, t, mu, eta)
    x1 = x(1);
    x2 = x(2);
    
    r = sqrt(x1^2 + x2^2);      %radial position  
                  
    gam = mu(1);
    gam1 = gam-1;

    gravity0 = mu(5);
    cp = mu(7);
    R = gam1*cp/gam;
    
%     rbot = mu(8);
%     pbot = mu(9);
%     Tbot = mu(10);
%     Ttop = mu(11);

    r0 = mu(12);
%     r1 = mu(13);

    m = mu(14);
    h = mu(15);
    c = mu(16);
    
%     latitude = 0;          %equatorial plane
%     declinationSun = 0;    %declination angle wrt the Sun
    
    %define local time
    cosChi = sqrt(0.5);
    absSinChi = sqrt(1-cosChi^2);
    
    %Computation F10.7 (let's assume it constant at first, the variation is at another scale)
    F10p7 = 100;
    F10p7_81 = 100;
    F10p7_mean = 0.5*(F10p7 + F10p7_81);
    
    rho = u(1);
    ru = u(2);
    rv = u(3);
    rE = u(4);
    
    r1 = 1.0/rho;
    uv = ru*r1;
    vv = rv*r1;
    ke = 0.5*(uv*uv+vv*vv);
    p = gam1*(rE-rho*ke);
    T = p/(R*rho);
    
    % Quantities
    gravity = gravity0*(r0^2/(r^2));
    H = R*T/(gravity);
    
    
    %Chapman integral 
    Rp = rho*H;
    Xp = r/H;
    y = sqrt(Xp/2)*abs(cosChi);
%     erfc = 1 - erf(y);
    
    Ierf = 0.5*(1+tanh(1000*(8-y)));
    a_erf = 1.06069630;
    b_erf = 0.55643831;
    c_erf = 1.06198960;
    d_erf = 1.72456090;
    f_erf = 0.56498823;
    g_erf = 0.06651874;

    erfcy = Ierf*(a_erf + b_erf*y)/(c_erf + d_erf*y + y*y) + (1-Ierf)*f_erf/(g_erf + y);
    
    
    IcosChi = 0.5*(1 + tanh(1000*cosChi));
    IsinChi = 0.5*(1 + tanh(1000*(r*absSinChi - r0)));
    
    alpha1 = Rp*erfcy*sqrt(0.5*pi*Xp);
    auxXp = (1-IcosChi)*IsinChi*Xp*(1-absSinChi);
    Rg = rho*H*exp(auxXp);
    alpha2 = (2*Rg - Rp*erfcy)*sqrt(0.5*pi*Xp);
    
    alpha = IcosChi*alpha1 + (1-IcosChi)*(IsinChi*alpha2 + (1-IsinChi)*1e32);
    
    crossSectionMean = 0;
    Q = 0;
    for iWave = 1:37
        lambda = eta(iWave);
        crossSection = eta(37+iWave);
        AFAC = eta(2*37+iWave);
        F74113 = eta(3*37+iWave);
        
%         tau = crossSection*alpha/m;
        
        slope0 = 1 + AFAC*(F10p7_mean-80);
        Islope = 0.5*(1+tanh(1000*(slope0-0.8)));
        slopeIntensity =  slope0*Islope + 0.8*(1-Islope);
        Intensity0 = F74113*slopeIntensity;
%         Intensity = Intensity0*exp(-tau);
        
        crossSectionMean = crossSectionMean + crossSection/37;
        Q = Q + crossSection*Intensity0/lambda;
    end
    
    Q = Q*exp(-crossSectionMean*alpha/m);
    
    eff0 = 0.6 - 5.54e-5*((r-r0)-65)^2;
    Ieff = 0.5*(1+tanh(1000*(eff0-0.2)));
    eff = 0.2 + (eff0-0.2)*Ieff;
    
    
    s_EUV = rho*eff*h*c*Q/m;