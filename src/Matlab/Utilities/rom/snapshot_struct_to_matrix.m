function [snapshots_matrix] = snapshot_struct_to_matrix(snapshots, pde, master)
    % Arrange snapshots struct into a matrix for POD
    n_train = length(snapshots);
    ncu = pde.ncu;
    ne  = pde.ne;
    npe = master.npe;
    
    snapshots_matrix = zeros(numel(snapshots{1}(:,1:ncu,:)), n_train);
    for i = 1:n_train
        tmp = snapshots{i}(:,1:ncu,:);
        snapshots_matrix(:,i) = reshape(tmp, [npe*ncu*ne, 1]);
    end
end