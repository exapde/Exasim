function HoverRT = enthalpyEval(speciesStruct, T, alpha)
    if nargin < 3; alpha = 1e8; end
    h_T0_to_T1 = nasa9eval_h(T, speciesStruct.a1, speciesStruct.b1);
    h_T1_to_T2 = nasa9eval_h(T, speciesStruct.a2, speciesStruct.b2);
    h_T2_to_Tmax = nasa9eval_h(T, speciesStruct.a3, speciesStruct.b3);
    HoverRT = switch1(alpha, speciesStruct.T1, T) * h_T0_to_T1 + switch2(alpha, speciesStruct.T1, speciesStruct.T2, T) * h_T1_to_T2 + switch3(alpha, speciesStruct.T2, T) * h_T2_to_Tmax;
end
