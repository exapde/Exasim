function lnkb_r = logBackwardCoefficients(Tb_r, lnkf_r, nu_f_kj, nu_b_kj, G0_j, P_atm, nr, ns)
    RU = 8.314471468617452;
    lnkb_r = zeros(nr,1,class(Tb_r));
    for ir = 1:nr %TODO: Tb_r should be a vector for 2T
        tmp = lnkf_r(ir); % TODO: Not correct for 2T models...need to evaluate with  Tbr
        tempTerm = RU * Tb_r;
        pressureTerm = log( P_atm / (tempTerm) );
        for k = 1:ns
            tmp = tmp + (nu_b_kj(k,ir) - nu_f_kj(k,ir)) * ( G0_j(k) - pressureTerm );
        end
        lnkb_r(ir) = tmp;
    end
end