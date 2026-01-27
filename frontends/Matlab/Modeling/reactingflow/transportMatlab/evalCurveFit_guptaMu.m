function out = evalCurveFit_guptaMu(T, params)
    % 0.1 exp( (A_n lnT + B_n) lnT + C_n )
    lnT = log(T);
    % return 0.1 * exp.((params.A * lnT .+ params.B) .* lnT + params.C)
    expTerm = exp(params.C);
    Tterm = T^(params.A * lnT + params.B);
    out = 0.1 * expTerm * Tterm;
end