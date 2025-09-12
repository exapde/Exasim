function out = switch3(alpha, T2, T)
    out = tanh( alpha * (T - T2) / pi ) * 0.5 + 0.5;
end