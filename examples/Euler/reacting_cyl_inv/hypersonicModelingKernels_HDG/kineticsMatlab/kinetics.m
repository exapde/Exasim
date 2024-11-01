function kinetics_params = kinetics()


% R_n = RU ./ Mw
% % rho_n = rho_inf * Y_vec
T_ref = 298.15;%mix.standardStateT()
P_atm = 101325.0;%mix.standardStateP()
A_r = [3.000e+16, 1.000e+16, 5.000e+09, 5.700e+06, 8.400e+06];
beta_r = [-1.60, -1.50, 0.00, 0.42, 0.00];
theta_r = [113200.0, 59360.0, 75500.0, 42938.0, 19400.0];
nu_f_kj = [
    0 0 0 1 0
    0 0 0 0 1
    0 0 1 0 0
    0 1 0 1 0
    0 1 1 0 0
];
nu_b_kj = [
    2 0 0 0 0
    0 2 0 0 0
    1 1 0 0 0
    1 0 1 0 0
    1 0 0 0 1
];
nu_f_kj = transpose(nu_f_kj); % NOTE: I manually messed up, went reaction by reaction, rather than species by species
nu_b_kj = transpose(nu_b_kj);
% alpha_jr = 
% [
%     1.0 1.0 0.23 0.23 0.23
%     1.0 1.0 0.20 0.20 0.20
%     22. 22. 22.  1.0  1.0
% ]
alpha_jr = [
    1.00 1.00 22.0
    1.00 1.00 22.0
    0.2333 0.20 22.0
    0.2333 0.20 1.00
    0.2333 0.20 1.00
];

ns = 5;%mix.nSpecies()
nr = 5;%mix.nReactions()

kinetics_params = struct("A_r", A_r, "beta_r", beta_r, "theta_r", theta_r, "nu_f_kj", nu_f_kj, "nu_b_kj", nu_b_kj, "alpha_jr", alpha_jr, "ns", ns, "nr", nr, "T_ref", T_ref, "P_atm", P_atm);

end
% gtst = [
%     -18.3451
%     -22.2584
%     -34.1779
%     -32.8287
%     -35.1012
% ]
% 



% mix.setState(rho_inf * Y_vec, 8000,1)

% %%%%%% Test
% thirdbody_r = zeros(5)
% kf_r = zeros(5); kb_r = zeros(5)
% lnkf_r = zeros(nr)
% lnkb_r = copy(lnkf_r)

% rhotst = [1.476033339116972e-03, 4.494900009825862e-04, 9.285632537926974e-08, 4.425397500274400e-06, 1.267967788517888e-09]
% retst = 7.417241701682229e+04

% % mix.setState(rhotst, retst,0)
% trange = 100:50:15000
% werr = zeros(length(trange), 5)
% for ti = 1:length(trange)
%     Tcurr = trange(ti)
%     mix.setState(Y_vec*rho_inf, Tcurr,1)
%     rho_i = mix.densities()
%     T = mix.T()
%     % rho_tilde = mix.densities() ./ mix.speciesMw()
%     % Gstd = G_n(mix.T())
%     % logForwardRateCoefficients!(lnkf_r, A_r, beta_r, theta_r, nr, mix.T())
%     % logBackwardCoefficients!(lnkb_r, mix.T(), lnkf_r, nu_f_kj, nu_b_kj, Gstd, RU, P_atm, nr)
%     % rateCoefficients!(kf_r,kb_r,lnkf_r,lnkb_r)
%     % thirdBodyFactor!(thirdbody_r, rho_tilde, alpha_jr, nr, ns)
%     % wtst, rr = netProductionRates(Mw, nu_f_kj, nu_b_kj, nr, ns, thirdbody_r, kf_r, kb_r, rho_tilde)
%     wtst, R_i = netProductionRatesTotal(rho_i, T, mix.speciesMw(), ns, nr)
%     werr(ti,:) = abs.((wtst - mix.netProductionRates()))
% end
% plot(trange, werr.+1e-13, yaxis=:log)

% %sym check
% @syms rho_sym[1:5], T_sym
% thirdbody_r_sym = zeros(Sym,5)
% kf_r_sym = copy(thirdbody_r_sym); 
% kb_r_sym = copy(thirdbody_r_sym); 
% lnkf_r_sym = copy(thirdbody_r_sym); 
% lnkb_r_sym = copy(thirdbody_r_sym)

% rho_tilde_sym = rho_sym ./ mix.speciesMw()
% Gsym = G_n(T_sym)
% logForwardRateCoefficients!(lnkf_r_sym, A_r, beta_r, theta_r, nr, T_sym)
% logBackwardCoefficients!(lnkb_r_sym, T_sym, lnkf_r_sym, nu_f_kj, nu_b_kj, Gsym, RU, P_atm, nr)
% rateCoefficients!(kf_r_sym, kb_r_sym, lnkf_r_sym, lnkb_r_sym)
% thirdBodyFactor!(thirdbody_r_sym, rho_tilde_sym, alpha_jr, nr, ns)
% wsym = netProductionRates(Mw, nu_f_kj, nu_b_kj, nr, ns, thirdbody_r_sym, kf_r_sym, kb_r_sym, rho_tilde_sym);

% function netProductionRatesTotal(rho_i, T, Mw, ns, nr)
% % TODO: could use a better name
%     Gformation = G_n(T)
%     thirdbody_r_sym = zeros(typeof(rho_i(1)), ns)
%     kf_r_sym = copy(thirdbody_r_sym); 
%     kb_r_sym = copy(thirdbody_r_sym); 
%     lnkf_r_sym = copy(thirdbody_r_sym); 
%     lnkb_r_sym = copy(thirdbody_r_sym)

%     rho_tilde = rho_i ./ Mw

%     logForwardRateCoefficients!(lnkf_r_sym, A_r, beta_r, theta_r, nr, T)
%     logBackwardCoefficients!(lnkb_r_sym, T, lnkf_r_sym, nu_f_kj, nu_b_kj, Gformation, RU, P_atm, nr)
%     rateCoefficients!(kf_r_sym, kb_r_sym, lnkf_r_sym, lnkb_r_sym)
%     thirdBodyFactor!(thirdbody_r_sym, rho_tilde, alpha_jr, nr, ns)
%     omega_i,R_i = netProductionRates(Mw, nu_f_kj, nu_b_kj, nr, ns, thirdbody_r_sym, kf_r_sym, kb_r_sym, rho_tilde)
%     return omega_i,R_i
% end
% println((wtst - mix.netProductionRates())./mix.netProductionRates())








% function logForwardRateCoefficients!(lnkf_r, A_r, beta_r, theta_r, nr, Tf_r)
%     for ir = 1:nr
%         lnkf_r(ir) = log(A_r(ir)) + beta_r(ir) * log(Tf_r(ir)) - theta_r(ir)/Tf_r(ir)
%     end
% end

% function logBackwardCoefficients!(lnkb_r, Tb_r, lnkf_r, nu_f_kj,nu_b_kj,G0_j, RU, P_atm, nr)
%     for ir = 1:nr
%         tmp = lnkf_r(ir) % TODO: NOT CORRECT..need to evaluate with  Tbr
%         tempTerm = RU * Tb_r(ir)
%         pressureTerm = log( P_atm / (tempTerm) )
%         for k = 1:ns
%             % println(G0_j(k) - pressureTerm )
%             tmp = tmp + (nu_b_kj(k,ir) - nu_f_kj(k,ir)) * ( G0_j(k) - pressureTerm )
%             % tmp = tmp + ( G0_j(k)/tempTerm - pressureTerm )
%         end
%         lnkb_r(ir) = tmp
%     end
% end

% function rateCoefficients!(kf_r,kb_r,lnkf_r,lnkb_r)
%     kf_r(:) = exp.(lnkf_r(:))
%     kb_r(:) = exp.(lnkb_r(:))
% end

% function thirdBodyFactor!(thirdbody_r, rho_tilde, alpha_jr, nr, ns)
%     % TODO: return 1 if not a third body
%     for ir = 1:nr
%         tmp = 0
%         for is = 1:ns
%             tmp = tmp + alpha_jr(is,ir) * rho_tilde(is)
%         end
%         TB_r(ir) = tmp
%     end
% end

% function netProductionRates!(omega_k, Mw, nu_f_kj, nu_b_kj, nr, ns, TB_r, kf_r, kb_r, rho_tilde)
%     % @variables R_r[1:nr] % Temporary storage of symbolic variables
%     % R_r = collect(R_r)
%     Rr = zeros(nr)
%     for r = 1:nr
%         Cf = 1
%         Cr = 1
%         for s = 1:ns
%             Cf = Cf * rho_tilde(s)^nu_f_kj(s,r)
%             Cr = Cr * rho_tilde(s)^nu_b_kj(s,r)
%         end
%         R_r(r) = (kf_r(r) * Cf - kb_r(r) * Cr) * TB_r(r)
%     end

%     for k = 1:ns
%         tmp = 0
%         for r = 1:nr
%             tmp = tmp + Mw(k) * (nu_b_kj(k,r) - nu_f_kj(k,r)) * R_r(r)
%         end
%         omega_k(k) = tmp
%     end
% end

