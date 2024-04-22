function tensorproduct(x::Array{Float64,2},porder::Int)::Tuple{Array{Float64,2},Array{Float64,2},Array{Float64,2},Array{Float64,2}}
# #TENSORPRODUCT Vandermonde matrices for tensor product polynomials in [0,1]^d
# #
# #   [F,FX,FY,FZ] = TENSORPRODUCT(X,P)
# #
# #      X:         Coordinates of the points wherethe polynomials
# #                 are to be evaluated (npoints,dim)
# #      PORDER:    Maximum order of the polynomials consider. That
# #                 is all polynomials of complete degree up to p,
# #                 npoly = (PORDER+1)*(PORDER+1)*(PORDER+1)
# #      F:         Vandermonde matrix (npoints,npoly)
# #      FX:        Vandermonde matrix for the derivative of the Koornwinder
# #                 polynomials w.r.t. x (npoints,npoly)
# #      FY:        Vandermonde matrix for the derivative of the Koornwinder
# #                 polynomials w.r.t. y (npoints,npoly)
# #      FZ:        Vandermonde matrix for the derivative of the Koornwinder
# #                 polynomials w.r.t. z (npoints,npoly)
# #

n,dim = size(x);


if dim== 1 # 1D
        f,fx = legendrepolynomial(x,porder);  # Legendre basis
        fy=[0.0 0.0];
        fz=[0.0 0.0];
elseif dim == 2 # 2D
        g1,gx = legendrepolynomial(x[:,1],porder); # Legendre basis in x direction
        g2,gy = legendrepolynomial(x[:,2],porder); # Legendre basis in y direction
        f  = zeros(Float64,n,(porder+1)*(porder+1));
        fx = 0*f;
        fy = 0*f;
        fz=[0.0 0.0];
        # perform tensor product to obtain the shape functions and their
        # derivatives on the unit square
        for ii=1:n
            tm = kron(g2[ii,:],g1[ii,:]);
            f[ii,:] =  kron(g2[ii,:],g1[ii,:]);
            fx[ii,:] = kron(g2[ii,:],gx[ii,:]);
            fy[ii,:] = kron(gy[ii,:],g1[ii,:]);
        end
elseif dim==3
        g1,gx=legendrepolynomial(x[:,1],porder); # Legendre basis in x direction
        g2,gy=legendrepolynomial(x[:,2],porder); # Legendre basis in y direction
        g3,gz=legendrepolynomial(x[:,3],porder); # Legendre basis in z direction
        f  = zeros(Float64,n,(porder+1)*(porder+1)*(porder+1));
        fx = 0*f;
        fy = 0*f;
        fz = 0*f;
        # perform tensor product to obtain the shape functions and their
        # derivatives on the unit cube
        for ii=1:n
            f[ii,:] =  kron(g3[ii,:],kron(g2[ii,:],g1[ii,:]));
            fx[ii,:] = kron(g3[ii,:],kron(g2[ii,:],gx[ii,:]));
            fy[ii,:] = kron(g3[ii,:],kron(gy[ii,:],g1[ii,:]));
            fz[ii,:] = kron(gz[ii,:],kron(g2[ii,:],g1[ii,:]));
        end
else
        error("Only can handle dim=1, dim=2 or dim=3");
end

return f, fx, fy, fz;

end
