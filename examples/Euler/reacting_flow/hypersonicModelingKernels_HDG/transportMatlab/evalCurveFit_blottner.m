function out = evalCurveFit_blottner(T, params)
    % 0.1 exp( (A_n lnT + B_n) lnT + C_n )
    lnT = log(T);
    out =  0.1 * exp.((params.A * lnT + params.B) .* lnT + params.C);
end
