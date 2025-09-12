
function [omega_k, R_r] = netProductionRates(Mw, nu_f_kj, nu_b_kj, nr, ns, thirdbody_r, kf_r, kb_r, rho_tilde)
    % R_r = collect(R_r)
    % R_r = zeros(typeof(kb_r(1)),nr)
    % omega_k = zeros(typeof(kb_r(1)), ns)
    R_r = zeros(nr,1,class(kb_r(1)));
    omega_k = zeros(ns,1,class(kb_r(1)));
    for r = 1:nr
        Cf = 1;
        Cr = 1;
        for s = 1:ns
            Cf = Cf * rho_tilde(s)^nu_f_kj(s,r);
            Cr = Cr * rho_tilde(s)^nu_b_kj(s,r);
        end
        % println(kb_r(r) * Cr)
        R_r(r) = (kf_r(r) * Cf - kb_r(r) * Cr) * thirdbody_r(r);
    end
    % println(R_r)

    for k = 1:ns
        tmp = 0;
        for r = 1:nr
            tmp = tmp + Mw(k) * (nu_b_kj(k,r) - nu_f_kj(k,r)) * R_r(r);
        end
        omega_k(k) = tmp;
    end
end