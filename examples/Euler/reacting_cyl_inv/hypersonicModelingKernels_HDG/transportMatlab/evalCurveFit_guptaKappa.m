function out = evalCurveFit_guptaKappa(T, params)
    % 0.1 exp( (A_n lnT + B_n) lnT + C_n )
    lnT = log(T);
    lnT2 = lnT*lnT;
    lnT3 = lnT2*lnT;
    % return 0.1 * exp.((params.A * lnT .+ params.B) .* lnT + params.C)
    expTerm = exp(params.E);
    Tterm = T^(params.A * lnT3 + params.B * lnT2 + params.C * lnT + params.D);
    out = expTerm * Tterm;
    out = 427.715996578 * out;
end