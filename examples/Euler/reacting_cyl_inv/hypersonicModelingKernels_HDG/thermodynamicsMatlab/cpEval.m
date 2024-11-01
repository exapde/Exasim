function CPoverR = cpEval(speciesStruct, T, alpha)
    if nargin < 3; alpha = 1e8; end
    cp_T0_to_T1 = nasa9eval_cp(T, speciesStruct.a1, speciesStruct.b1);
    cp_T1_to_T2 = nasa9eval_cp(T, speciesStruct.a2, speciesStruct.b2);
    cp_T2_to_Tmax = nasa9eval_cp(T, speciesStruct.a3, speciesStruct.b3);
    CPoverR = switch1(alpha, speciesStruct.T1, T) * cp_T0_to_T1 + switch2(alpha, speciesStruct.T1, speciesStruct.T2, T) * cp_T1_to_T2 + switch3(alpha, speciesStruct.T2, T) * cp_T2_to_Tmax;
end
