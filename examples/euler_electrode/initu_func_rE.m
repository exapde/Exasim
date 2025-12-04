function u0 = initu_func_rE(dgnodes, param)
    z_tilde = dgnodes(:,1,:);
    r_tilde = dgnodes(:,2,:);

    u0 = 2.5 + 2.5*(3-1)/(1 + exp(3*(r_tilde-2)));
    % u0=u0*10;
end
