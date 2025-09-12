function out = evalCurveFit_guptaDij(T, params)
    % exp(D) T^(A (lnT)^2 + B lnT + C)
    lnT = log(T);
    expD = exp(params.D);
    Tterm = T^(params.A * lnT^2 + params.B * lnT + params.C);
    out = expD * Tterm;
end