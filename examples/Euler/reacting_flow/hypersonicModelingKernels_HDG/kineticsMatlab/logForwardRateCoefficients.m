function lnkf_r = logForwardRateCoefficients(A_r, beta_r, theta_r, nr, Tf_r)
    lnkf_r = zeros(nr,1,class(Tf_r));
    for ir = 1:nr %TODO: Tf_r should be a vector for 2T
        lnkf_r(ir) = log(A_r(ir)) + beta_r(ir) * log(Tf_r) - theta_r(ir)/Tf_r;
    end
end