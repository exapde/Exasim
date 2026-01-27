function eos = equationofstate(T, rhos, rhoe)

    [species_thermo_structs, Mw, RU] = thermodynamicsModels;
    ns = length(Mw);
    rho_tilde = rhos ./ Mw;

    alpha = 1e4;
    fT = elementaryfunctions(T);
    
    eos = -sum(rho_tilde);
    for i = 1:ns
        Tsw1 = species_thermo_structs{i}.T1;
        Tsw2 = species_thermo_structs{i}.T2;
        fsw = switchfunctions(T, Tsw1, Tsw2, alpha);      
        c1 = nasa9_hcoeff(species_thermo_structs{i}.a1, species_thermo_structs{i}.b1);
        c2 = nasa9_hcoeff(species_thermo_structs{i}.a2, species_thermo_structs{i}.b2);
        c3 = nasa9_hcoeff(species_thermo_structs{i}.a3, species_thermo_structs{i}.b3);
        h1 = sum(c1.*fT);
        h2 = sum(c2.*fT);
        h3 = sum(c3.*fT);  
        H = fsw(1)*h1 + fsw(2)*h2 + fsw(3)*h3;      
        eos = eos + rho_tilde(i) * H;
    end
    eos = eos * T;
    eos = eos - rhoe / RU;

end

