function s_EUV = EUVsource(u, x, t, mu, eta)
    x1 = x(1);
    x2 = x(2);
    
    gam = mu(1);
    gam1 = gam - 1.0;
    Minf = mu(4);
    M2 = Minf^2;
    
    Fr2 = mu(5);
    omega = mu(6);
    Q0 = mu(12);
    M0 = mu(13);
    
    R0 = mu(10);
    R1 = mu(11);
    Ldim = mu(14);
    
    %% computation  
    z = sqrt(x1^2 + x2^2);      %radial position  
              
    sgnx1 = tanh(10000*x1);
    sgnx2 = tanh(10000*x2);
    sgnx2pow2 = tanh(10000*x2*x2);
    sgnx1x2 = tanh(10000*x1*x2);
    longitude = pi - pi/2*(1+sgnx1)*(1-sgnx2pow2) - pi/4*(2+sgnx1)*sgnx2 - sgnx1x2*atan((abs(x1) - abs(x2))/(abs(x1) + abs(x2)));
    
%     latitude = 0;          %equatorial plane
%     declinationSun = 0;    %declination angle wrt the Sun
    
    %define local time
    long_offset = omega*t;
    localTime = longitude + long_offset;
%     cosChi = sin(declinationSun)*sin(latitude) + cos(declinationSun)*cos(latitude)*cos(localTime);
    cosChi = cos(localTime);
    
    absSinChi = sqrt(1-cosChi^2);
    
    %Computation F10.7 (let's assume it constant at first, the variation is at another scale)
    F10p7 = 100;
    F10p7_81 = 100;
    F10p7_mean = 0.5*(F10p7 + F10p7_81);
    
    r = u(1);
    rho = exp(r);
    T = u(4)/sqrt(rho); 

    % Quantities
    gravity = Fr2*(R0^2/(z^2));
    H = T/(gam*M2*gravity);
    
    %Chapman integral 
    Rp = rho*H;
    Xp = z/H;
    y = sqrt(Xp/2)*abs(cosChi);
    
    Ierf = 0.5*(1+tanh(1000*(8-y)));
    a_erf = 1.06069630;
    b_erf = 0.55643831;
    c_erf = 1.06198960;
    d_erf = 1.72456090;
    f_erf = 0.56498823;
    g_erf = 0.06651874;

    erfcy = Ierf*(a_erf + b_erf*y)/(c_erf + d_erf*y + y*y) + (1-Ierf)*f_erf/(g_erf + y);
    
%     erfy = erf(y);
%     erfcy = exp(y^2)*(1-erfy);
    
    IcosChi = 0.5*(1 + tanh(100000*cosChi));
    IsinChi = 0.5*(1 + tanh(100000*(z*absSinChi - R0)));
    
    alpha1 = Rp*erfcy*sqrt(0.5*pi*Xp);
    auxXp = (1-IcosChi)*IsinChi*Xp*(1-absSinChi);
    Rg = rho*H*exp(auxXp);
    alpha2 = (2*Rg - Rp*erfcy)*sqrt(0.5*pi*Xp);
    
    alpha = IcosChi*alpha1 + (1-IcosChi)*(IsinChi*alpha2 + (1-IsinChi)*1e32);
    
    Q = 0;
    for iWave = 1:37
        lambda = eta(iWave);
        crossSection = eta(37+iWave);
        AFAC = eta(2*37+iWave);
        F74113 = eta(3*37+iWave);
        
        tau = M0*crossSection*alpha;
        
        slope0 = 1 + AFAC*(F10p7_mean-80);
        Islope = 0.5*(1+tanh(1000*(slope0-0.8)));
        slopeIntensity =  slope0*Islope + 0.8*(1-Islope);
        Intensity0 = F74113*slopeIntensity;
        Intensity = Intensity0*exp(-tau);
        
        Q = Q + crossSection*Intensity/lambda;
    end
    
%     H0 = 65/500*(R1-R0);
%     eff0 = 0.6 - 5.54e-11*(Ldim*((z-R0)-H0)).^2;
%     Ieff = 0.5*(1+tanh(1000*(eff0-0.2)));
%     eff = 0.2 + (eff0-0.2).*Ieff;
    eff = mu(15);
    
    s_EUV = gam*gam1*M2*eff*Q0*Q;