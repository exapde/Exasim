
function omega_i = netProductionRatesTotal(rho_i, T, Mw, kinetics_params, species_thermo_structs)
    % TODO: could use a better name
    Gformation = G_n(T, species_thermo_structs);
    A_r = kinetics_params.A_r;
    beta_r = kinetics_params.beta_r;
    theta_r = kinetics_params.theta_r;
    nu_f_kj = kinetics_params.nu_f_kj;
    nu_b_kj = kinetics_params.nu_b_kj;
    alpha_jr = kinetics_params.alpha_jr;
    nr = kinetics_params.nr;
    ns = kinetics_params.ns;
    T_ref = kinetics_params.T_ref;
    P_atm = kinetics_params.P_atm;

    % thirdbody_r_sym = zeros(typeof(rho_i(1)), ns)
    % kf_r_sym = copy(thirdbody_r_sym); 
    % kb_r_sym = copy(thirdbody_r_sym); 
    % lnkf_r_sym = copy(thirdbody_r_sym); 
    % lnkb_r_sym = copy(thirdbody_r_sym)

    rho_tilde = rho_i ./ Mw;

    lnkf_r_sym = logForwardRateCoefficients(A_r, beta_r, theta_r, nr, T);
    lnkb_r_sym = logBackwardCoefficients(T, lnkf_r_sym, nu_f_kj, nu_b_kj, Gformation, P_atm, nr, ns);
    [kf_r_sym, kb_r_sym] = rateCoefficients(lnkf_r_sym, lnkb_r_sym);
    thirdbody_r_sym = thirdBodyFactor(rho_tilde, alpha_jr, nr, ns);
    [omega_i,~] = netProductionRates(Mw, nu_f_kj, nu_b_kj, nr, ns, thirdbody_r_sym, kf_r_sym, kb_r_sym, rho_tilde);
end
