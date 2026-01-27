function thirdbody_r = thirdBodyFactor(rho_tilde, alpha_jr, nr, ns)
    thirdbody_r = zeros(nr,1,class(rho_tilde(1)));
    for ir = 1:3
        tmp = 0;
        for is = 1:ns
            tmp = tmp + alpha_jr(is,ir) * rho_tilde(is);
        end
        thirdbody_r(ir) = tmp;
    end
    for ir = 4:5
        thirdbody_r(ir) = 1.0;
    end
end
