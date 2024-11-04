function out = switch1(alpha, T1, T)
    out = tanh( -alpha * (T - T1) / pi ) * 0.5 + 0.5;
end
