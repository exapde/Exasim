function pgf = cmeshmap3(pf, dgnodes, lw, ll)

xr = dgnodes(:,1,:); 
yr = dgnodes(:,2,:); 
[xgf,ygf] = cmeshmap2( pf(:,1), pf(:,2), xr(:), yr(:), lw, ll);

pgf = 0*dgnodes;
pgf(:,1,:) = reshape(xgf, size(xr));
pgf(:,2,:) = reshape(ygf, size(xr));

