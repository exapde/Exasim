function G = nasa9eval_G(T, a, b)
    T2 = T * T;
    T3 = T2 * T;
    T4 = T3 * T;
    Tinv = 1.0 / T;
    logT = log(T);
    % println(-a(1) * 1.0/(2*T2))
    G = (-a(1) * 1.0/(2*T2) + a(2) * (logT + 1.0) * Tinv + a(3) * (1 - logT) - a(4) * T/2 - a(5) * T2 / 6.0 - a(6) * T3 / 12.0 - a(7) * T4 / 20.0 + b(1)/T - b(2));
end

