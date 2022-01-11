#ifndef __MATVEC
#define __MATVEC

#include "ioutilities.cpp"
void MatVec(dstype *w, solstruct &sol, resstruct &res, appstruct &app, masterstruct &master,
      meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, dstype *v, 
      dstype *u, dstype *Ru, Int backend)
{   
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)    
    Int npe = common.npe; // number of nodes on master element    
    Int ne = common.ne1; // number of elements in this subdomain 
    Int nd = common.nd;
    Int N = npe*ncu*ne;
    Int order = common.matvecOrder;
    dstype epsilon = common.matvecTol;
#ifdef HAVE_ENZYME
//TODO: there might be a cleaner way to do this...matvecAD v. matvecFD functions? 
    ArrayInsert(sol.dudg, v, npe, nc, ne, 0, npe, 0, ncu, 0, ne, backend); //insert v into dudgs
    dResidual(sol, res, app, master, mesh, tmp, common, handle, backend);
    ArrayAXPBY(w, u, res.dRu, 0.0, 1, N, backend);
#else
    if (order==1) {
        // calculate w = u + epsilon*v
        ArrayAXPBY(w, u, v, 1.0, epsilon, N, backend);
                            
        // insert (u+epsilon*v) into udg
        ArrayInsert(sol.udg, w, npe, nc, ne, 0, npe, 0, ncu, 0, ne, backend);  

        // compute the residual R(u+epsilon*v)
        Residual(sol, res, app, master, mesh, tmp, common, handle, backend);

        // calculate w = J(u)*v = (R(u+epsilon*v)-R(u))/epsilon    
        ArrayAXPBY(w, res.Ru, Ru, 1.0/epsilon, -1.0/epsilon, N, backend);
    }
    else if (order==2) {
        // calculate w = u - epsilon*v
        ArrayAXPBY(w, u, v, 1.0, -epsilon, N, backend);

        // insert (u-epsilon*v) into udg
        ArrayInsert(sol.udg, w, npe, nc, ne, 0, npe, 0, ncu, 0, ne, backend);  

        // compute the residual R(u-epsilon*v)
        Residual(sol, res, app, master, mesh, tmp, common, handle, backend);

        // copy res.Ru to Ru
        ArrayCopy(Ru, res.Ru, N, backend);
        
        // calculate w = u + epsilon*v
        ArrayAXPBY(w, u, v, 1.0, epsilon, N, backend);

        // insert (u+epsilon*v) into udg
        ArrayInsert(sol.udg, w, npe, nc, ne, 0, npe, 0, ncu, 0, ne, backend);  

        // compute the residual R(u+epsilon*v)
        Residual(sol, res, app, master, mesh, tmp, common, handle, backend);
        
        // calculate w = J(u)*v = (R(u+epsilon*v)-R(u-epsilon*v))/(2*epsilon)    
        ArrayAXPBY(w, res.Ru, Ru, 0.5/epsilon, -0.5/epsilon, N, backend);
    }
    else
        error("Matrix-vector multiplication order is not implemented");
#endif
}
#endif


