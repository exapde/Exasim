/*
    uresidual.cpp

    This file contains functions for computing element and face residuals in a LDG discretization framework, 
    including support for time-dependent problems, source terms, and automatic differentiation via Enzyme.

    Functions:

    - RuElemBlock: Computes the element residuals for a block of elements, including source and flux contributions.
    - RuElem: Loops over element blocks and calls RuElemBlock for each block.
    - dRuElemBlock (HAVE_ENZYME): Computes the derivative of element residuals using Enzyme for automatic differentiation.
    - dRuElem (HAVE_ENZYME): Loops over element blocks and calls dRuElemBlock for each block.

    - RuFaceBlock: Computes the face residuals for a block of faces, handling both interior and boundary faces.
    - RuFace: Loops over face blocks and calls RuFaceBlock for each block.
    - dRuFaceBlock (HAVE_ENZYME): Computes the derivative of face residuals using Enzyme for automatic differentiation.
    - dRuFace (HAVE_ENZYME): Loops over face blocks and calls dRuFaceBlock for each block.

    Key Concepts:
    - Element and face residuals are computed at Gauss points and projected to nodes.
    - Supports time-dependent problems, source terms, and custom flux/source functions via driver routines.
    - Uses temporary arrays for intermediate computations and supports CUDA/cuBLAS for backend acceleration.
    - Optional debug output via EXADEBUG macro.
    - Enzyme support for automatic differentiation (HAVE_ENZYME macro).

    Arguments (common across functions):
    - solstruct &sol: Solution data structure.
    - resstruct &res: Residual data structure.
    - appstruct &app: Application-specific data.
    - masterstruct &master: Master element/face data.
    - meshstruct &mesh: Mesh connectivity and geometry.
    - tempstruct &tmp: Temporary storage for computations.
    - commonstruct &common: Common parameters and settings.
    - cublasHandle_t handle: cuBLAS handle for GPU computations.
    - Int e1, e2, f1, f2: Element/face range indices.
    - Int backend: Backend identifier (CPU/GPU).
    - Int ib: Face type (interior/boundary).

    Note:
    - Functions with HAVE_ENZYME macro are only compiled when Enzyme AD is enabled.
    - Debug output is controlled by the EXADEBUG macro.
*/
#ifndef __URESIDUAL
#define __URESIDUAL

void RuElemBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
        Int e1, Int e2, Int backend)
{        
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg) 
    Int ncs = common.ncs;// number of compoments of (sdg) 
    Int ncw = common.ncw;// number of compoments of (wdg) 
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe; // number of nodes on master element
    Int nge = common.nge; // number of gauss points on master element        

    Int ne = e2-e1;
    Int nn =  npe*ne; 
    Int nga = nge*ne;   
    //Int n0 = 0;                        // xg
    Int n1 = nga*ncx;                  // Xx
    Int n2 = nga*(ncx+nd*nd);          // jac        
    Int n5 = 0;                        // fg    
    Int n4 = nga*ncu*nd;               // sg  
    Int n3 = nga*(ncu*nd+ncu);         // ug
    Int n6 = nga*(ncu*nd+ncu+nc);      // wg
    Int nm = nge*e1*(ncx+nd*nd+1);
    
    // udg = tmp.tempg[n3] at gauss points on element
    //GetElemNodes(tmp.tempn, sol.udg, npe, nc, 0, nc, e1, e2, backend);
    GetArrayAtIndex(tmp.tempn, sol.udg, &mesh.eindudg1[npe*nc*e1], nn*nc);
    Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, master.shapegt, nge, npe, ne*nc, backend);
    if (ncw>0) {
        GetElemNodes(tmp.tempn, sol.wdg, npe, ncw, 0, ncw, e1, e2);    
        Node2Gauss(handle, &tmp.tempg[n6], tmp.tempn, master.shapegt, nge, npe, ne*ncw, backend);
    }
    
    if (common.tdep) { // for time-dependent problem                
        // calculate sdg = sdg-udg*dtfactor
        ArrayAXPBY(&tmp.tempg[n4], &sol.sdgg[nge*ncs*e1], &tmp.tempg[n3], one, -common.dtfactor, nga*ncu);            
        
        if (common.tdfunc==1) {
            // calculate the time derivative function Tdfunc(xdg, udg, odg)
            TdfuncDriver(&tmp.tempg[n5], &sol.elemg[nm], &tmp.tempg[n3], &sol.odgg[nge*nco*e1], 
                &tmp.tempg[n6], mesh, master, app, sol, tmp, common, nge, e1, e2, backend);

            // calculate (sdg-udg*dtfactor)*Tdfunc(xdg, udg, odg) 
            ArrayAXY(&tmp.tempg[n4], &tmp.tempg[n4], &tmp.tempg[n5], one, nga*ncu);                
        }
        
        if (common.source==1) {            
            // calculate the source term Source(xdg, udg, odg, wdg)
            SourceDriver(&tmp.tempg[n5], &sol.elemg[nm], &tmp.tempg[n3], &sol.odgg[nge*nco*e1], 
                &tmp.tempg[n6], mesh, master, app, sol, tmp, common, nge, e1, e2, backend);

            // calculate Source(xdg, udg, odg) + (sdg-udg*dtfactor)*Tdfunc(xdg, udg, odg) 
            ArrayAXPBY(&tmp.tempg[n4], &tmp.tempg[n4], &tmp.tempg[n5], one, one, nga*ncu);            
        }        
    }
    else {        
        // calculate the source term Source(xdg, udg, odg, wdg)
        SourceDriver(&tmp.tempg[n4], &sol.elemg[nm], &tmp.tempg[n3], &sol.odgg[nge*nco*e1],  
            &tmp.tempg[n6], mesh, master, app, sol, tmp, common, nge, e1, e2, backend);                 
    }                
        
    FluxDriver(&tmp.tempg[n5], &sol.elemg[nm], &tmp.tempg[n3], &sol.odgg[nge*nco*e1],  
            &tmp.tempg[n6], mesh, master, app, sol, tmp, common, nge, e1, e2, backend);    

    // RuSource(&res.Rue[npe*ncu*e1], &tmp.tempg[n4], &sol.elemg[nm+n2], master.shapegw, nge, npe, ncu, ne);
    // RuFlux(&res.Rue[npe*ncu*e1], &tmp.tempg[n5], &sol.elemg[nm+n1], &master.shapegw[npe*nge], nge, npe, ncu, nd, ne);

    // calculate sum_j Flux_j(u) * Xx(:,j,i)  at gauss points on element
    //         rg = sg, fg, Xx, jac
    ApplyXx4(&tmp.tempg[n3], &tmp.tempg[n4], &tmp.tempg[n5], &sol.elemg[nm+n1], &sol.elemg[nm+n2],
            nge, nd, ncu, ne);
    
    //(Source(u) + timesource, w)_K + (Flux(u), nabla w)_K = sum_j (Flux_j(u), dw/ dx_j)_K  = sum_j (jac Flux_j(u), dw/dxi_i dxi_i/dx_j)_T
    // = dw/dxi_i * (sum_j Flux_j(u) * Xx(:,j,i)) 
    Gauss2Node(handle, &res.Rue[npe*ncu*e1], &tmp.tempg[n3], master.shapegw, nge*(nd+1), npe, ncu*ne, backend);                
      
#ifdef EXADEBUG                       
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

#ifdef HAVE_ENZYME
//// Method 2
void dRuElemBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
        Int e1, Int e2, Int backend)
{        
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg) 
    Int ncs = common.ncs;// number of compoments of (sdg) 
    Int ncw = common.ncw;// number of compoments of (wdg) 
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe; // number of nodes on master element
    Int nge = common.nge; // number of gauss points on master element        

    Int ne = e2-e1;
    Int nn =  npe*ne; 
    Int nga = nge*ne;   
    Int n1 = nga*ncx;                  // Xx
    Int n2 = nga*(ncx+nd*nd);          // jac        
    Int n5 = 0;                        // fg    
    Int n4 = nga*ncu*nd;               // sg  
    Int n3 = nga*(ncu*nd+ncu);         // ug
    Int n6 = nga*(ncu*nd+ncu+nc);      // wg
    Int n7 = nga*(ncu*nd+ncu+nc+ncw);               // dfg
    Int n8 = nga*(ncu*nd+ncu+nc+ncw+ncu*nd);        // dsg  
    Int n9 = nga*(ncu*nd+ncu+nc+ncw+ncu*nd+ncu);    // dug
    Int n0 = nga*(ncu*nd+ncu+nc+ncw+ncu*nd+ncu+nc); // dwg
    Int nm = nge*e1*(ncx+nd*nd+1);

    GetArrayAtIndex(tmp.tempn, sol.udg, &mesh.eindudg1[npe*nc*e1], nn*nc);
    Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, master.shapegt, nge, npe, ne*nc, backend);        
    GetArrayAtIndex(tmp.tempn, sol.dudg, &mesh.eindudg1[npe*nc*e1], nn*nc);
    Node2Gauss(handle, &tmp.tempg[n9], tmp.tempn, master.shapegt, nge, npe, ne*nc, backend);        


    ArraySetValue(&tmp.tempg[n7], 0.0, nga*ncu*nd); // so df = 0.
    ArraySetValue(&tmp.tempg[n8], 0.0, nga*ncu);
    if (ncw>0) {
        GetElemNodes(tmp.tempn, sol.wdg, npe, ncw, 0, ncw, e1, e2);    
        Node2Gauss(handle, &tmp.tempg[n6], tmp.tempn, master.shapegt, nge, npe, ne*ncw, backend);
        GetElemNodes(tmp.tempn, sol.dwdg, npe, ncw, 0, ncw, e1, e2);    
        Node2Gauss(handle, &tmp.tempg[n0], tmp.tempn, master.shapegt, nge, npe, ne*ncw, backend);
    }

    if (common.tdep) { // for time-dependent problem                
        ArrayAXPBY(&tmp.tempg[n4], &sol.sdgg[nge*ncs*e1], &tmp.tempg[n3], one, -common.dtfactor, nga*ncu);

        // ArrayCopy(&tmp.tempg[n8], &tmp.tempg[n9], nga*ncu, backend);
        // ArrayMultiplyScalar(&tmp.tempg[n8], -common.dtfactor, nga*ncu, backend); 
        ArrayAXPBY(&tmp.tempg[n8], &tmp.tempg[n8], &tmp.tempg[n9], one, -common.dtfactor, nga*ncu);
        
        if (common.tdfunc==1) {
            // calculate the time derivative function Tdfunc(xdg, udg, odg)
            TdfuncDriver(&tmp.tempg[n5], &sol.elemg[nm], &tmp.tempg[n3], &sol.odgg[nge*nco*e1], 
                &tmp.tempg[n6], mesh, master, app, sol, tmp, common, nge, e1, e2, backend);

            // calculate (sdg-udg*dtfactor)*Tdfunc(xdg, udg, odg) 
            ArrayAXY(&tmp.tempg[n4], &tmp.tempg[n4], &tmp.tempg[n5], one, nga*ncu);                
            // calculate dsdg*Tdfunc(xdg, udg, odg) 
            ArrayAXY(&tmp.tempg[n8], &tmp.tempg[n8], &tmp.tempg[n5], one, nga*ncu);                
        }
        
        if (common.source==1) {            
            // calculate the source term Source(xdg, udg, odg, wdg) and dSource
            SourceDriver(&tmp.tempg[n5], &tmp.tempg[n7], &sol.elemg[nm], &tmp.tempg[n3], &tmp.tempg[n9], &sol.odgg[nge*nco*e1], 
                &tmp.tempg[n6], &tmp.tempg[n0], mesh, master, app, sol, tmp, common, nge, e1, e2, backend);
            
            // calculate Source(xdg, udg, odg) + (sdg-udg*dtfactor)*Tdfunc(xdg, udg, odg) 
            ArrayAXPBY(&tmp.tempg[n4], &tmp.tempg[n4], &tmp.tempg[n5], one, one, nga*ncu);            
            // calculate dSource(xdg, udg, odg) + dsdg*Tdfunc(xdg, udg, odg) 
            ArrayAXPBY(&tmp.tempg[n8], &tmp.tempg[n8], &tmp.tempg[n7], one, one, nga*ncu);            
        }        
    }
    else {        
        // calculate the source term Source(xdg, udg, odg, wdg)
        SourceDriver(&tmp.tempg[n4], &tmp.tempg[n8], &sol.elemg[nm], &tmp.tempg[n3], &tmp.tempg[n9], &sol.odgg[nge*nco*e1],  
            &tmp.tempg[n6], &tmp.tempg[n0], mesh, master, app, sol, tmp, common, nge, e1, e2, backend);             
    }    


    FluxDriver(&tmp.tempg[n5], &tmp.tempg[n7], &sol.elemg[nm], &tmp.tempg[n3], &tmp.tempg[n9], &sol.odgg[nge*nco*e1], &sol.dodgg[nge*nco*e1],
            &tmp.tempg[n6], &tmp.tempg[n0], mesh, master, app, sol, tmp, common, nge, e1, e2, backend);   

    // calculate sum_j dFlux_j(u) * Xx(:,j,i)  at gauss points on element
    // .       rg = sg, fg, Xx, jac
    ApplyXx4(&tmp.tempg[n9], &tmp.tempg[n8], &tmp.tempg[n7], &sol.elemg[nm+n1], &sol.elemg[nm+n2],
        nge, nd, ncu, ne);
    Gauss2Node(handle, &res.dRue[npe*ncu*e1], &tmp.tempg[n9], master.shapegw, nge*(nd+1), npe, ncu*ne, backend);             

#ifdef EXADEBUG                       
    writearray2file(common.fileout + "EnzymeRuElem_uge.bin", &tmp.tempg[n3], nge*nc*ne, backend);  
    writearray2file(common.fileout + "EnzymeRuElem_fge.bin", &tmp.tempg[n4], nge*ncu*nd*ne, backend);  
    writearray2file(common.fileout + "EnzymeRuElem_rne.bin", tmp.tempn, npe*ncu*ne, backend);
    writearray2file(common.fileout + "EnzymeRuElem_rqe.bin", res.Rue, npe*ncu*common.ne1, backend);
    writearray2file(common.fileout + "EnzymedRuElem_rqe.bin", res.dRue, npe*ncu*common.ne1, backend);
#endif                  
}

void dRuElem(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common,         
        cublasHandle_t handle, Int nbe1, Int nbe2, Int backend)
{    
    for (Int j=nbe1; j<nbe2; j++) {
        Int e1 = common.eblks[3*j]-1;
        Int e2 = common.eblks[3*j+1];            
        dRuElemBlock(sol, res, app, master, mesh, tmp, common, handle, e1, e2, backend);
    }                     
}
#endif                  


// Calculate Ruf = <fhat(xdg, uhat, udg, odg, nl), w>_F 
void RuFaceBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, 
        cublasHandle_t handle, Int f1, Int f2, Int ib, Int backend)
{            
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int ncw = common.ncw;
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
    Int n4 = nga*(ncu);                         // ug1
    Int n5 = nga*(ncu+nc);                      // wg1
    Int n6 = nga*(ncu+nc+ncw);                  // ug2
    Int n7 = nga*(ncu+nc+ncw+nc);               // wg2
    Int n8 = nga*(ncu+2*nc+2*ncw);              // fhg
    //Int n7 = nga*(ncu+2*nc+ncw);                // wdg
    Int nm = ngf*f1*(ncx+nd+1);
    
    // uhg = tmp.tempg[n3] at gauss points on face
    GetElemNodes(tmp.tempn, sol.uh, npf, ncu, 0, ncu, f1, f2);
    
    // ug1 = tmp.tempg[n4] at gauss points on face
    //GetFaceNodes(tmp.tempn, sol.udg, mesh.facecon, npf, nc, npe, nc, f1, f2, 1, backend);      
    GetArrayAtIndex(&tmp.tempn[nn*ncu], sol.udg, &mesh.findudg1[npf*nc*f1], nn*nc);    
    if (ncw>0)
        GetFaceNodes(&tmp.tempn[nn*(ncu+nc)], sol.wdg, mesh.facecon, npf, ncw, npe, ncw, f1, f2, 1);      
    
    if (ib==0) {
        // ug2 = tmp.tempg[n6] at gauss points on face
        //GetFaceNodes(tmp.tempn, sol.udg, mesh.facecon, npf, nc, npe, nc, f1, f2, 2, backend);     
        GetArrayAtIndex(&tmp.tempn[nn*(ncu+nc+ncw)], sol.udg, &mesh.findudg2[npf*nc*f1], nn*nc);
        if (ncw>0)
            GetFaceNodes(&tmp.tempn[nn*(ncu+nc+ncw+nc)], sol.wdg, mesh.facecon, npf, ncw, npe, ncw, f1, f2, 2);        
        Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, master.shapfgt, ngf, npf, nf*(ncu+2*nc+2*ncw), backend);
    }
    else
        Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, master.shapfgt, ngf, npf, nf*(ncu+nc+ncw), backend);
        
    // calculate fhat
    if (ib==0) { // interior faces                
        FhatDriver(&tmp.tempg[n8], &sol.faceg[nm+n0], &tmp.tempg[n4], &tmp.tempg[n6], 
        &sol.og1[ngf*nco*f1], &sol.og2[ngf*nco*f1], &tmp.tempg[n5], &tmp.tempg[n7], 
        &tmp.tempg[n3], &sol.faceg[nm+n1], mesh, master, app, sol, tmp, common, ngf, f1, f2, backend);      
    }
    else { // boundary faces      
        FbouDriver(&tmp.tempg[n8], &sol.faceg[nm+n0], &tmp.tempg[n4], &sol.og1[ngf*nco*f1], 
                &tmp.tempg[n5], &tmp.tempg[n3], &sol.faceg[nm+n1], mesh, master, app, 
                sol, tmp, common, ngf, f1, f2, ib, backend);        
    }        
            
    // evaluate fhg * jac at gauss points on face
    ApplyJacFhat(&tmp.tempg[n3], &tmp.tempg[n8], &sol.faceg[nm+n2], nga, ncu, ngf);    
    
    // <fhat, w>_F = <jac fhat, w>_T = w * (fhg * jac): npf*ncu*nf        
    Gauss2Node(handle, &res.Rh[npf*ncu*f1], &tmp.tempg[n3], master.shapfgw, ngf, npf, nf*ncu, backend);            
    
#ifdef EXADEBUG                           
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

#ifdef HAVE_ENZYME
//// Method 2
// Calculate Ruf = <fhat(xdg, uhat, udg, odg, nl), w>_F 
void dRuFaceBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, 
        cublasHandle_t handle, Int f1, Int f2, Int ib, Int backend)
{            
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int ncw = common.ncw;
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
    Int n4 = nga*(ncu);                         // ug1
    Int n5 = nga*(ncu+nc);                      // wg1
    Int n6 = nga*(ncu+nc+ncw);                  // ug2
    Int n7 = nga*(ncu+nc+ncw+nc);               // wg2
    Int n8 = nga*(ncu+2*nc+2*ncw);              // fhg
    Int n9 = nga*(ncu+2*nc+2*ncw+2*ncu*nd);       //duhg
    Int n10 = nga*(ncu+2*nc+2*ncw+2*ncu*nd+ncu);        // dug1
    Int n11 = nga*(ncu+2*nc+2*ncw+2*ncu*nd+ncu+nc);     // dwg1
    Int n12 = nga*(ncu+2*nc+2*ncw+2*ncu*nd+ncu+nc+ncw); // dug2
    Int n13 = nga*(ncu+2*nc+2*ncw+2*ncu*nd+ncu+2*nc+ncw); // dwg2
    Int n14 = nga*(ncu+2*nc+2*ncw+2*ncu*nd+ncu+2*nc+2*ncw); // dfhg

    Int nm = ngf*f1*(ncx+nd+1);

    Int nnd = nn*(ncu+2*nc+2*ncw);

    ////// First, u, uL, uR on faces are loaded 
    // uhg = tmp.tempg[n3] at gauss points on face
    GetElemNodes(tmp.tempn, sol.uh, npf, ncu, 0, ncu, f1, f2, backend);

    // ug1 = tmp.tempg[n4] at gauss points on face
    GetArrayAtIndex(&tmp.tempn[nn*ncu], sol.udg, &mesh.findudg1[npf*nc*f1], nn*nc, backend);    
    if (ncw>0)
        GetFaceNodes(&tmp.tempn[nn*(ncu+nc)], sol.wdg, mesh.facecon, npf, ncw, npe, ncw, f1, f2, 1, backend);      
    
    if (ib==0) {
        // ug2 = tmp.tempg[n6] at gauss points on face
        GetArrayAtIndex(&tmp.tempn[nn*(ncu+nc+ncw)], sol.udg, &mesh.findudg2[npf*nc*f1], nn*nc, backend);
        if (ncw>0)
            GetFaceNodes(&tmp.tempn[nn*(ncu+nc+ncw+nc)], sol.wdg, mesh.facecon, npf, ncw, npe, ncw, f1, f2, 2, backend);        
        Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, master.shapfgt, ngf, npf, nf*(ncu+2*nc+2*ncw), backend);
    }
    else
        Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, master.shapfgt, ngf, npf, nf*(ncu+nc+ncw), backend);

    ///// Repeat above steps to get du, duL, duR; indices for face information should be the same
    // TODO: tmp.tmpn was used to store node information above. Should I overwrite parts of tmp.tmpn
    //       or should I append new information on to the end of it?
    //
    //       For example, should the next line of code be tmp.tempn or tmp.tempn[nn*(ncu+2*nc+2*ncw)]  

    // duh = tmp.tempg[n9] at gauss points on face
    GetElemNodes(&tmp.tempn[nnd], sol.duh, npf, ncu, 0, ncu, f1, f2, backend);

    // dug1 = tmp.tempg[n10] at gauss points on face
    //GetFaceNodes(tmp.tempn, sol.udg, mesh.facecon, npf, nc, npe, nc, f1, f2, 1, backend);      
    GetArrayAtIndex(&tmp.tempn[nnd+nn*ncu], sol.dudg, &mesh.findudg1[npf*nc*f1], nn*nc, backend);    
    if (ncw>0)
        GetFaceNodes(&tmp.tempn[nnd+nn*(ncu+nc)], sol.dwdg, mesh.facecon, npf, ncw, npe, ncw, f1, f2, 1, backend);      
    
    if (ib==0) {
        // dug2 = tmp.tempg[n12] at gauss points on face
        //GetFaceNodes(tmp.tempn, sol.udg, mesh.facecon, npf, nc, npe, nc, f1, f2, 2, backend);     
        GetArrayAtIndex(&tmp.tempn[nnd+nn*(ncu+nc+ncw)], sol.dudg, &mesh.findudg2[npf*nc*f1], nn*nc, backend);
        if (ncw>0)
            GetFaceNodes(&tmp.tempn[nnd+nn*(ncu+nc+ncw+nc)], sol.dwdg, mesh.facecon, npf, ncw, npe, ncw, f1, f2, 2, backend);        
        Node2Gauss(handle, &tmp.tempg[n9], &tmp.tempn[nnd], master.shapfgt, ngf, npf, nf*(ncu+2*nc+2*ncw), backend);
    }
    else
        Node2Gauss(handle, &tmp.tempg[n9], &tmp.tempn[nnd], master.shapfgt, ngf, npf, nf*(ncu+nc+ncw), backend);

//// Apply Enzyme on fhat or fbou
    // calculate fhat
    if (ib==0) { // interior faces                
        FhatDriver(&tmp.tempg[n8], &tmp.tempg[n14], &sol.faceg[nm+n0], &tmp.tempg[n4], &tmp.tempg[n10], &tmp.tempg[n6], &tmp.tempg[n12],
        &sol.og1[ngf*nco*f1], &sol.dog1[ngf*nco*f1], &sol.og2[ngf*nco*f1], &sol.dog2[ngf*nco*f1], &tmp.tempg[n5], &tmp.tempg[n11], &tmp.tempg[n7], &tmp.tempg[n13],
        &tmp.tempg[n3], &tmp.tempg[n9], &sol.faceg[nm+n1], mesh, master, app, sol, tmp, common, ngf, f1, f2, backend);      
    }
    else { // boundary faces      
        FbouDriver(&tmp.tempg[n8], &tmp.tempg[n14], &sol.faceg[nm+n0], &tmp.tempg[n4], &tmp.tempg[n10], &sol.og1[ngf*nco*f1], 
                &sol.dog1[ngf*nco*f1], &tmp.tempg[n5], &tmp.tempg[n11], &tmp.tempg[n3], &tmp.tempg[n9], &sol.faceg[nm+n1], mesh, master, app, 
                sol, tmp, common, ngf, f1, f2, ib, backend);        
    } 
    // evaluate dfhg * jac
    ApplyJac(&tmp.tempg[n9], &tmp.tempg[n14], &sol.faceg[nm+n2], nga, ncu, ngf, backend);
    Gauss2Node(handle, &res.dRh[npf*ncu*f1], &tmp.tempg[n9], master.shapfgw, ngf, npf, nf*ncu, backend);   
#ifdef EXADEBUG                           
    writearray2file(common.fileout + NumberToString(ib) + "RuFace_uhgf.bin", &tmp.tempg[n3], ngf*ncu*nf, backend);  
    writearray2file(common.fileout + NumberToString(ib) + "RuFace_fgf.bin", &tmp.tempg[n4], ngf*ncu*nf, backend);  
    writearray2file(common.fileout + NumberToString(ib) + "RuFace_rnf.bin", tmp.tempn, npf*ncu*nf, backend);
    writearray2file(common.fileout + NumberToString(ib) + "RuFace_ruf.bin", res.Ruf, npe*ncu*common.ne1, backend);
#endif              
}

void dRuFace(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common,
        cublasHandle_t handle, Int nbf1, Int nbf2, Int backend)
{    
    for (Int j=nbf1; j<nbf2; j++) {
        Int f1 = common.fblks[3*j]-1;
        Int f2 = common.fblks[3*j+1];    
        Int ib = common.fblks[3*j+2];    
        dRuFaceBlock(sol, res, app, master, mesh, tmp, common, handle, f1, f2, ib, backend);
    }                          
}
#endif

#endif

