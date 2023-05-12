function [Wmodes, svals] = snapshot_POD(Phi, W, pde, master)
%Performs POD on snapshot matrix and returns POD modes and 
    if nargin < 2
        % Scalar weight
        W = 1;
        pde = [];
        master = [];
    end
    WPhi = 0*Phi;
    n_train = size(Phi, 2);
    for i = 1:n_train
        WPhi(:,i) = apply_weight_to_vector(W, Phi(:,i), pde, master);
    end
    [V,D] = eig(Phi'*Phi);
    [~,ii] = sort(diag(D),'descend');
    D = D(ii,ii);
    V = V(:,ii);
    Wmodes = Phi*V;
    svals = diag(D);
end