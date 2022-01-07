#ifndef __URESIDUAL
#define __URESIDUAL

#include "../application/fluxDriver.cpp"
#include "../application/sourceDriver.cpp"
#include "../application/tdfuncDriver.cpp"
#include "../application/fhatDriver.cpp"
#include "../application/fbouDriver.cpp"
#include "../application/avfdDriver.cpp"

void RuElemBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
        Int e1, Int e2, Int backend)
{        
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg) 
    Int ncs = common.ncs;// number of compoments of (sdg) 
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe; // number of nodes on master element
    Int nge = common.nge; // number of gauss points on master element        

    Int ne = e2-e1;
    Int nn =  npe*ne; 
    Int nga = nge*ne;   
    Int n0 = 0;                                 // xg
    Int n1 = nga*ncx;                           // Xx
    Int n2 = nga*(ncx+nd*nd);                   // jac        
    Int n5 = 0;                        // fg    
    Int n4 = nga*ncu*nd;               // sg  
    Int n3 = nga*(ncu+ncu*nd);         // ug
    Int nm = nge*e1*(ncx+nd*nd+1);
    
    // udg = tmp.tempg[n3] at gauss points on element
    //GetElemNodes(tmp.tempn, sol.udg, npe, nc, 0, nc, e1, e2, backend);
    GetArrayAtIndex(tmp.tempn, sol.udg, &mesh.eindudg1[npe*nc*e1], nn*nc, backend);
    Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, master.shapegt, nge, npe, ne*nc, backend);
        
    if (common.tdep) { // for time-dependent problem                
        // calculate sdg = sdg-udg*dtfactor
        ArrayAXPBY(&tmp.tempg[n4], &sol.sdgg[nge*ncs*e1], &tmp.tempg[n3], one, -common.dtfactor, nga*ncu, backend);            
        
        if (common.tdfunc==1) {
            // calculate the time derivative function Tdfunc(xdg, udg, odg)
            TdfuncDriver(&tmp.tempg[n5], &sol.elemg[nm], &tmp.tempg[n3], &sol.odgg[nge*nco*e1], 
                mesh, master, app, sol, tmp, common, nge, e1, e2, backend);

            // calculate (sdg-udg*dtfactor)*Tdfunc(xdg, udg, odg) 
            ArrayAXY(&tmp.tempg[n4], &tmp.tempg[n4], &tmp.tempg[n5], one, nga*ncu, backend);                
        }
        
        if (common.source==1) {
            // calculate the source term Source(xdg, udg, odg)
            SourceDriver(&tmp.tempg[n5], &sol.elemg[nm], &tmp.tempg[n3], &sol.odgg[nge*nco*e1], 
                mesh, master, app, sol, tmp, common, nge, e1, e2, backend);

            // calculate Source(xdg, udg, odg) + (sdg-udg*dtfactor)*Tdfunc(xdg, udg, odg) 
            ArrayAXPBY(&tmp.tempg[n4], &tmp.tempg[n4], &tmp.tempg[n5], one, one, nga*ncu, backend);            
        }        
    }
    else {        
        // calculate the source term Source(xdg, udg, odg)
        SourceDriver(&tmp.tempg[n4], &sol.elemg[nm], &tmp.tempg[n3], &sol.odgg[nge*nco*e1],  
            mesh, master, app, sol, tmp, common, nge, e1, e2, backend);        
    }    
        
    FluxDriver(&tmp.tempg[n5], &sol.elemg[nm],  &tmp.tempg[n3], &sol.odgg[nge*nco*e1],  
            mesh, master, app, sol, tmp, common, nge, e1, e2, backend);    
    
    // calculate sum_j Flux_j(u) * Xx(:,j,i)  at gauss points on element
    // .       rg = sg, fg, Xx, jac
    ApplyXx4(&tmp.tempg[n3], &tmp.tempg[n4], &tmp.tempg[n5], &sol.elemg[nm+n1], &sol.elemg[nm+n2],
            nge, nd, ncu, ne, backend);
    
    //(Source(u) + timesource, w)_K + (Flux(u), nabla w)_K = sum_j (Flux_j(u), dw/ dx_j)_K  = sum_j (jac Flux_j(u), dw/dxi_i dxi_i/dx_j)_T
    // = dw/dxi_i * (sum_j Flux_j(u) * Xx(:,j,i)) 
    Gauss2Node(handle, &res.Rue[npe*ncu*e1], &tmp.tempg[n3], master.shapegw, nge*(nd+1), npe, ncu*ne, backend);                
        
#ifdef DEBUG                       
    writearray2file(common.fileout + "RuElem_uge.bin", &tmp.tempg[n3], nge*nc*ne, backend);  
    writearray2file(common.fileout + "RuElem_fge.bin", &tmp.tempg[n4], nge*ncu*nd*ne, backend);  
    writearray2file(common.fileout + "RuElem_rne.bin", tmp.tempn, npe*ncu*ne, backend);
    writearray2file(common.fileout + "RuElem_rqe.bin", res.Rue, npe*ncu*common.ne1, backend);
#endif                  
}

void RuElem(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common,         
        cublasHandle_t handle, Int nbe1, Int nbe2, Int backend)
{    
    for (Int j=nbe1; j<nbe2; j++) {
        Int e1 = common.eblks[3*j]-1;
        Int e2 = common.eblks[3*j+1];    
        RuElemBlock(sol, res, app, master, mesh, tmp, common, handle, e1, e2, backend);
    }                     
}

// Calculate Ruf = <fhat(xdg, uhat, udg, odg, nl), w>_F 
void RuFaceBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, 
        cublasHandle_t handle, Int f1, Int f2, Int ib, Int backend)
{            
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe; // number of nodes on master element
    Int npf = common.npf; // number of nodes on master face           
    Int ngf = common.ngf; // number of gauss poInts on master face              

    Int nf = f2-f1;
    Int nn =  npf*nf; 
    Int nga = ngf*nf;   
    Int n0 = 0;                                 // xg
    Int n1 = nga*ncx;                           // nlg
    Int n2 = nga*(ncx+nd);                      // jac
    Int n3 = nga*(0);                           // uhg    
    Int n5 = nga*(ncu);                         // ug1
    Int n6 = nga*(ncu+nc);                      // ug2
    Int n4 = nga*(ncu+2*nc);                    // fhg
    Int nm = ngf*f1*(ncx+nd+1);
    
    // uhg = tmp.tempg[n3] at gauss points on face
    GetElemNodes(tmp.tempn, sol.uh, npf, ncu, 0, ncu, f1, f2, backend);
    
    // ug1 = tmp.tempg[n5] at gauss points on face
    //GetFaceNodes(tmp.tempn, sol.udg, mesh.facecon, npf, nc, npe, nc, f1, f2, 1, backend);      
    GetArrayAtIndex(&tmp.tempn[nn*ncu], sol.udg, &mesh.findudg1[npf*nc*f1], nn*nc, backend);    
    
    if (ib==0) {
        // ug2 = tmp.tempg[n6] at gauss points on face
        //GetFaceNodes(tmp.tempn, sol.udg, mesh.facecon, npf, nc, npe, nc, f1, f2, 2, backend);     
        GetArrayAtIndex(&tmp.tempn[nn*(ncu+nc)], sol.udg, &mesh.findudg2[npf*nc*f1], nn*nc, backend);        
        Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, master.shapfgt, ngf, npf, nf*(ncu+2*nc), backend);
    }
    else
        Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, master.shapfgt, ngf, npf, nf*(ncu+nc), backend);
    
    // calculate fhat
    if (ib==0) { // interior faces                
        FhatDriver(&tmp.tempg[n4], &sol.faceg[nm+n0], &tmp.tempg[n3], &tmp.tempg[n5], 
        &tmp.tempg[n6], &sol.og1[ngf*nco*f1], &sol.og2[ngf*nco*f1], &sol.faceg[nm+n1], 
        mesh, master, app, sol, tmp, common, ngf, f1, f2, backend);      
    }
    else { // boundary faces      
        FbouDriver(&tmp.tempg[n4], &sol.faceg[nm+n0], &tmp.tempg[n3], &tmp.tempg[n5], 
        &sol.og1[ngf*nco*f1], &sol.faceg[nm+n1], mesh, master, app, sol, tmp, common, 
        ngf, f1, f2, ib, backend);
    }        
        
    // evaluate fhg * jac at gauss points on face
    ApplyJac(&tmp.tempg[n3], &tmp.tempg[n4], &sol.faceg[nm+n2], nga, ncu, ngf, backend);    
    
    // <fhat, w>_F = <jac fhat, w>_T = w * (fhg * jac): npf*ncu*nf        
    Gauss2Node(handle, &res.Rh[npf*ncu*f1], &tmp.tempg[n3], master.shapfgw, ngf, npf, nf*ncu, backend);            
    
#ifdef DEBUG                           
    writearray2file(common.fileout + NumberToString(ib) + "RuFace_uhgf.bin", &tmp.tempg[n3], ngf*ncu*nf, backend);  
    writearray2file(common.fileout + NumberToString(ib) + "RuFace_fgf.bin", &tmp.tempg[n4], ngf*ncu*nf, backend);  
    writearray2file(common.fileout + NumberToString(ib) + "RuFace_rnf.bin", tmp.tempn, npf*ncu*nf, backend);
    writearray2file(common.fileout + NumberToString(ib) + "RuFace_ruf.bin", res.Ruf, npe*ncu*common.ne1, backend);
#endif              
}

void RuFace(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common,
        cublasHandle_t handle, Int nbf1, Int nbf2, Int backend)
{    
    for (Int j=nbf1; j<nbf2; j++) {
        Int f1 = common.fblks[3*j]-1;
        Int f2 = common.fblks[3*j+1];    
        Int ib = common.fblks[3*j+2];    
        RuFaceBlock(sol, res, app, master, mesh, tmp, common, handle, f1, f2, ib, backend);
    }                          
}

#endif

