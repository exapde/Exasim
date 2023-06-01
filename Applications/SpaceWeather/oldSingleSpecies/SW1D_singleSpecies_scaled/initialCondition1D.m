function UDG = initialCondition1D(mesh,UDG,mu)
    
    Fr = mu(4);

    Tbot = mu(8);
    Ttop = mu(9);
    R0 = mu(10);
    Ldim = mu(14);
    h0 = 35000/Ldim;

    a0 = (-1 + Fr^2*R0);

% Radial position
nElems = size(mesh.dgnodes,3);
nNodesElem = size(mesh.dgnodes,1);


for iElem = 1:nElems
    for iNode = 1:nNodesElem
        r = mesh.dgnodes(iNode,1,iElem);

        T = Ttop - (Ttop-Tbot)*exp(-(r-R0)/h0);
        logp_p0 = a0*h0/Ttop*log(1+Ttop/Tbot*(exp((r-R0)/h0)-1));
        rtilde = logp_p0 - log(T);
        rho = exp(rtilde);
        srT = sqrt(rho)*T;
        
        UDG(iNode,1,iElem) = rtilde;
        UDG(iNode,3,iElem) = srT;
    end
end