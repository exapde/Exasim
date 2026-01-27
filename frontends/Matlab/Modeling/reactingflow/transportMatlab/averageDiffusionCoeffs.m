function D_n = averageDiffusionCoeffs(T, X_i, Y_i, Mw, p, gupta_structs)
    ns = length(X_i);
    D_n = zeros(ns,1, class(X_i(1)));

    Dij = fDij(T, gupta_structs) / p * 10.1325; %dimensionalization
    
    for i = 1:ns
        denom = 0.0 * X_i(1);
        for j = 1:ns
            if i ~= j
                denom = denom + X_i(j) / Dij(i,j);
            end
        end
        D_n(i) = (1.0 - X_i(i)) / denom;
    end
end