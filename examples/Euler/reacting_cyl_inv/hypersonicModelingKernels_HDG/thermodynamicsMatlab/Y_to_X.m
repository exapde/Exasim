function out = Y_to_X(Y_i, Mw)
    X = Y_i ./ Mw;
    X = X ./ sum(X);
    out =  X;
end
