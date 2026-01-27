function h = nasa9eval_h(T, a, b)
    T2 = T * T;
    T3 = T2 * T;
    T4 = T3 * T;
    Tinv = 1.0 / T;
    logT = log(T);
    h = -a(1) * 1/T2  + a(2) * logT * Tinv  + a(3) + a(4) * T/2  + a(5) * T2 / 3 + a(6) * T3 / 4 + a(7) * T4/5 + b(1) / T ;
end
