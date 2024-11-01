function out = switch2(alpha, T1, T2, T)
    out = 1.0 - switch1(alpha, T1, T) - switch3(alpha, T2, T);
end
