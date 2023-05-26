function [out] = subsample_vector(V, pde, master,ii)
    npe = master.npe;
    ncu = pde.ncu;
    ne = pde.ne;
    u_shape = [npe ncu ne];
    u_flat = [npe*ncu*ne 1];
    V = reshape(V, u_shape);
    out_tmp = V(:,:,ii);
    out = reshape(out_tmp, [npe*ncu*length(ii) 1]);
end