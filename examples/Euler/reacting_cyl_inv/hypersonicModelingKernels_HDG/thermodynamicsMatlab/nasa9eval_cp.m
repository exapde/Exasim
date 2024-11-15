function cp = nasa9eval_cp(T, a, b)
    T2 = T * T;
    T3 = T2 * T;
    T4 = T3 * T;
    Tinv = 1.0 / T;
    logT = log(T);
    cp = a(1) * 1/T2 + a(2) * Tinv + a(3) + a(4) * T + a(5) * T2 + a(6) * T3  + a(7) * T4;
end