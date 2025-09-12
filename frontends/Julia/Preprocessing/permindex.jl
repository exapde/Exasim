function permindex(plocfc,dim,elemtype)

npf = size(plocfc,1);

if dim==1
    ind = 1;
elseif dim==2
    ind = collect(npf:-1:1);
elseif dim==3
    if elemtype==0
        ind = zeros(Int,npf,3);

        # [1 3 2]
        plocfc2 = copy(plocfc);
        plocfc2[:,1] = plocfc[:,2];
        plocfc2[:,2] = plocfc[:,1];
        ind[:,1] = xiny(round.(plocfc,digits=8),round.(plocfc2,digits=8),1);

        # [2 1 3]
        plocfc2 = copy(plocfc);
        plocfc2[:,1] = 1.0 .- plocfc[:,1] .- plocfc[:,2];
        ind[:,2] = xiny(round.(plocfc,digits=8),round.(plocfc2,digits=8),1);

        # [3 2 1]
        plocfc2 = copy(plocfc);
        plocfc2[:,2] = 1.0 .- plocfc[:,1] .- plocfc[:,2];
        ind[:,3] = xiny(round.(plocfc,digits=8),round.(plocfc2,digits=8),1);
    else
        ind = zeros(Int,npf,4);

        # [1 4 3 2]
        plocfc2 = copy(plocfc);
        plocfc2[:,1] = plocfc[:,2];
        plocfc2[:,2] = plocfc[:,1];
        ind[:,1] = xiny(round.(plocfc,digits=8),round.(plocfc2,digits=8),1);

        # [2 1 4 3]
        plocfc2 = copy(plocfc);
        plocfc2[:,2] = 1.0 .- plocfc[:,2];
        ind[:,2] = xiny(round.(plocfc,digits=8),round.(plocfc2,digits=8),1);

        # [3 2 1 4]
        plocfc2 = copy(plocfc);
        plocfc2[:,1] = 1.0 .- plocfc[:,2];
        plocfc2[:,2] = 1.0 .- plocfc[:,1];
        ind[:,3] = xiny(round.(plocfc,digits=8),round.(plocfc2,digits=8),1);

        # [4 3 2 1]
        plocfc2 = copy(plocfc);
        plocfc2[:,1] = 1.0 .- plocfc[:,1];
        ind[:,4] = xiny(round.(plocfc,digits=8),round.(plocfc2,digits=8),1);
    end
end

return ind

end
