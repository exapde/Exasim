function [u_rb, u_rb_history] = run_steady_ROM(mu, u_0, Phi, pde, mesh, master, dmd, opt,ii,JPhi_train)
% System to solve is W^T R(W u_rb) = 0 or || R(W u_rb) ||
% W = reduced basis
% u_rb = generalized coordinates

% TODOs:  
% - weight with Cholesky factorization of matrix; simpler to keep track of
% - probably make separate functions for LSPG, Galerkin
% - need appropriate measures of convergence for each 
% - make a rom data structure to hold all parameters

    pde.physicsparam = mu;
    fileapp = "datain/app.bin";
    pde.runmode = 20; % run mode > 5: R(u,mu) and J(u,mu) v
    writeapp(pde,fileapp,'native');

    maxiter = 20;
    nltol = 1e-8;
    iter = 1;
    galerkin_flag = 0; lspg_flag = 0; 

    if opt == "galerkin" || opt == "Galerkin" || opt == "g" || opt=="G"
        galerkin_flag = 1;
    end
    if opt == "lspg" || opt == "LSPG" || opt == "L" || opt == "lspg"
        lspg_flag = 1;
    end


    npe = master.npe;
    ncu = pde.ncu;
    ne = pde.ne;
    u_shape = [npe ncu ne];
    u_flat = [npe*ncu*ne 1];
    p_rb = size(Phi, 2);

    u_rb = u_0;
    Phiu = Phi * u_rb;
    mesh.udg = reshape(Phiu, [master.npe, pde.ncu, pde.ne]);
    mesh.dudg = 0*mesh.udg;
    writesol(pde,mesh,master,dmd)
    runcode(pde);

    R0 = reshape(getsolution('dataout/out_Ru_test',dmd,npe),u_flat);
    JPhi = 0*R0;

    if galerkin_flag
        b = -Phi' * R0;
    elseif lspg_flag 
        b = R0' * R0; %TODO: should be weighted norm
        JPhi = 0*R0;
        WJPhi = 0*JPhi;
        [Mi, M] = massinv(pde, master,mesh);
        weight = 1;
    end

    %TODO: flag for recording history
    %   include u_rb, tolerance, etc. 
    u_rb_history = cell(maxiter,1);
    u_rb_history{1} = u_rb;

    [ne, p_rb_2, n_train] = size(JPhi_train);
    ne_hr = length(ii);
    JPhi_train_subsampled = zeros(ne_hr*ncu*npe, size(JPhi_train, 2), size(JPhi_train, 3));
    JPhi_subsampled = zeros(ne_hr*ncu*npe, p_rb);
    for j = 1:p_rb_2
        for k = 1:n_train
            JPhi_train_subsampled(:,j,k) = subsample_vector(JPhi_train(:,j,k), pde, master, ii);
        end
    end


    while norm(b) > nltol && iter < maxiter
        Phiu = Phi * u_rb;
        mesh.udg = reshape(Phiu, u_shape);
        for i = 1:p_rb
            % Evaluate JPhi by calling Exasim's Jv against each column of Phi
            mesh.dudg = reshape(Phi(:,i), u_shape);
            writesol(pde,mesh,master,dmd)
            runcode(pde);
            Jv = reshape(getsolution('dataout/out_Jv_test',dmd,npe),u_flat);
            JPhi(:,i) = Jv;
        end
%         JPhi_subsampled = JPhi(ii,:);
%         JPhi_subsampled = subsample_vector(JPhi, pde, master, ii);
        for j = 1:p_rb_2
            JPhi_subsampled(:,j) = subsample_vector(JPhi(:,j), pde, master, ii);
        end
%         JPhi_train_subsampled = JPhi_train(ii,:,:);
%         for k = 1:n_train
%             JPhi_train_subsampled = 
%         end
%         [ne_hr, p_rb_2, n_train] = size(JPhi_train_subsampled);

        A_hr = reshape(JPhi_train_subsampled, [ne_hr*ncu*npe*p_rb_2 n_train]);
        b_hr = JPhi_subsampled(:);

        theta = (A_hr'*A_hr)\(A_hr'*b_hr);
        JPhi = 0;
        for k = 1:n_train
            JPhi = JPhi + theta(k) * JPhi_train(:,:,k);
        end

        R = reshape(getsolution('dataout/out_Ru_test',dmd,npe),u_flat);
    %%%%%%%%% Galerkin
    if galerkin_flag
        A =  Phi' * JPhi;
        b = -Phi' * R;
    end
    if lspg_flag
%         for j = 1:p_rb
%             JPhi_subsampled(:,j) = subsample_vector(JPhi(:,j), pde, master, ii);
%             Phi_subsampled(:,j) = subsample_vector(Phi(:,j), pde, master, ii);
%         end
% %         Phi_subsampled = subsample_vector(Phi, pde, master, ii);
%         R_subsampled = subsample_vector(R, pde, master, ii);
%         A =  Phi(ii,:)' * JPhi(ii,:);
%         b = -Phi(ii,:)' * R(ii,:);
        A =  JPhi' * JPhi;
        b = -JPhi' * R;

%         A =  JPhi' * WJPhi;
%         b = -(WJPhi)' * R;
    end
    %%%%%%%%%
    
        du_rb = A\b;
        u_rb = u_rb + du_rb;
        u_rb_history{iter+1} = u_rb;
        iter = iter+1;
        disp(norm(b));
    end
end