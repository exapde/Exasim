function sca = chemreactingeval(u,str,rho_scale,u_scale,T)
%EULERVAL Calculates derived quantities for the Euler equation variables.
%   SCA=EULERVAL(U,STR,GAM)
%
%      UP(npl,4,nt):   np plus states
%      STR:            String used to specify requested quantity
%                      - STR: 'r' Density
%                      - STR: 'u' u_x velocity
%                      - STR: 'v' u_y velocity
%                      - STR: 'p' Density
%                      - STR: 'M' Density
%                      - STR; 's' Entropy
%      GAM:            Value of Gamma
%      SCA(npl,4,nt):  Scalar field requested by STR 
%

[thermo_structs, Mw, ~] = thermodynamicsModels();
ns = length(Mw);

rho_s = u(:,1:ns,:);
rho_i = rho_s * rho_scale;

if strcmp(str,'r')
    sca = sum(rho_i,2);
elseif strcmp(str,'p')
    sca = pressure(T, rho_i, Mw);
elseif strcmp(str,'gam')
    X = X_i(rho_i, Mw);
    Y = Y_i(rho_i);
    sca = mixtureFrozenGamma(T, Mw, Y, X, thermo_structs);
elseif strcmp(str,'c')
    X = X_i(rho_i, Mw);
    Y = Y_i(rho_i);
    gam = mixtureFrozenGamma(T, Mw, Y, X, thermo_structs);
    p = pressure(T, rho_i, Mw);
    density = sum(rho_i,2);
    sca = sqrt(gam .* p ./ density);
elseif strcmp(str,'M')
    X = X_i(rho_i, Mw);
    Y = Y_i(rho_i);
    gam = mixtureFrozenGamma(T, Mw, Y, X, thermo_structs);
    p = pressure(T, rho_i, Mw);
    density = sum(rho_i,2);
    c = sqrt(gam .* p ./ density);
    rhou = u(:,ns+1,:) * (rho_scale * u_scale);
    rhov = u(:,ns+2,:) * (rho_scale * u_scale);
    uv = rhou./density;
    vv = rhov./density;
    u2 = sqrt(uv.^2+vv.^2);
    sca = u2./c;
elseif strcmp(str,'s')
    X = X_i(rho_i, Mw);
    Y = Y_i(rho_i);
    gam = mixtureFrozenGamma(T, Mw, Y, X, thermo_structs);
    p = pressure(T, rho_i, Mw);
    density = sum(rho_i,2);
    sca = p./(density.^gam);
elseif strcmp(str,'u')
    rhou = u(:,ns+1,:) * (rho_scale * u_scale);
    density = sum(rho_i,2);
    sca = rhou./density;
elseif strcmp(str,'v')
    rhov = u(:,ns+2,:) * (rho_scale * u_scale);
    density = sum(rho_i,2);
    sca = rhov./density;
elseif strcmp(str,'h')
    p = pressure(T, rho_i, Mw);
    density = sum(rho_i,2);
    sca = u(:,ns+3,:)./density+p./density;    
else
    error('unknonw case');
end

end

function p = pressure(T, r_i, Mw)
    RU = 8.314471468617452;
    p = T .* sum(r_i ./ (Mw(:)'),2) * RU;
end

function out = X_i(rho_i, Mw)
    conc = rho_i./(Mw(:)');
    out =  conc ./ sum(conc,2);
end

function out = Y_i(rho_i)
    out =  rho_i ./ sum(rho_i,2);
end

function cp = nasa9eval_cp(T, a, b)
    T2 = T .* T;
    T3 = T2 .* T;
    T4 = T3 .* T;
    Tinv = 1.0 ./ T;
    cp = a(1) * 1./T2 + a(2) * Tinv + a(3) + a(4) * T + a(5) * T2 + a(6) * T3  + a(7) * T4;
end

% function h = nasa9eval_h(T, a, b)
%     T2 = T .* T;
%     T3 = T2 .* T;
%     T4 = T3 .* T;
%     Tinv = 1.0 ./ T;
%     logT = log(T);
%     h = (-a(1) * 1./T2  + a(2) * logT .* Tinv  + a(3) + a(4) * T/2  + a(5) * T2 / 3 + a(6) * T3 / 4 + a(7) * T4/5 + b(1) ./ T );
% end

function CPoverR = cpEval(speciesStruct, T, alpha)
    if nargin < 3; alpha = 1e8; end
    cp_T0_to_T1 = nasa9eval_cp(T, speciesStruct.a1, speciesStruct.b1);
    cp_T1_to_T2 = nasa9eval_cp(T, speciesStruct.a2, speciesStruct.b2);
    cp_T2_to_Tmax = nasa9eval_cp(T, speciesStruct.a3, speciesStruct.b3);
    CPoverR = switch1(alpha, speciesStruct.T1, T) .* cp_T0_to_T1 + switch2(alpha, speciesStruct.T1, speciesStruct.T2, T) .* cp_T1_to_T2 + switch3(alpha, speciesStruct.T2, T) .* cp_T2_to_Tmax;
end

function out = cp_n(T, species_thermo_structs)
    npe = size(T,1);
    ne = size(T,3);
    ns = length(species_thermo_structs);
    out = zeros(npe, ns, ne);
    for i = 1:ns
        out(:,i,:) = cpEval(species_thermo_structs{i}, T);
    end
end

function cp_m = getCpsMass(T, Mw, species_thermo_structs)
    RU = 8.314471468617452;
    cp = cp_n(T, species_thermo_structs);
    cp_m = cp * RU ./ (Mw(:)');
end

function out = mixtureFrozenCpMass(T, Mw, Y, species_thermo_structs)
    cp = getCpsMass(T, Mw, species_thermo_structs);
    out = sum(cp .* Y, 2);
end

function cv_m = getCvsMass(T, Mw, species_thermo_structs)
    RU = 8.314471468617452;
    cp = cp_n(T, species_thermo_structs);
    cv_m = (cp - 1.0) * RU ./ (Mw(:)');
end

function out = mixtureFrozenCvMass(T, Mw, Y, thermo_structs)
    cv = getCvsMass(T, Mw, thermo_structs);
    out =  sum(cv .* Y, 2);
end

function out = mixtureMw(Mw, X)
    out = sum(X .* (Mw(:)'), 2);
end

function out = mixtureFrozenCpMole(T, Mw, Y, X, thermo_structs)
    out =  mixtureFrozenCpMass(T, Mw, Y, thermo_structs) .* mixtureMw(Mw, X);
end

function out = mixtureFrozenCvMole(T, Mw, Y, X, thermo_structs)
    out =  mixtureFrozenCvMass(T, Mw, Y, thermo_structs) .* mixtureMw(Mw, X);
end

function out = mixtureFrozenGamma(T, Mw, Y, X, thermo_structs)
    cp = mixtureFrozenCpMole(T, Mw, Y, X, thermo_structs);
    cv = mixtureFrozenCvMole(T, Mw, Y, X, thermo_structs);
    out = cp ./ cv;
end

