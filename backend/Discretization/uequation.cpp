/*
  uequation.cpp

  This file contains functions for assembling and solving the HDG (Hybridizable Discontinuous Galerkin) system for a general PDE discretization. The main routines are responsible for computing element and face contributions, assembling global matrices, and handling Schur complement reduction for efficient solution of the HDG system.

  Functions:

  - uEquationElemBlock:
    Assembles element-wise contributions to the HDG system, including source terms, fluxes, and their derivatives. Handles time-dependent and steady-state problems, and supports auxiliary variables (e.g., w-equation coupling).

  - uEquationElemFaceBlock:
    Assembles face-wise contributions, including numerical fluxes and boundary/interface conditions. Handles the mapping between element and face degrees of freedom, and supports auxiliary variables and coupled boundary conditions.

  - uEquationSchurBlock:
    Performs Schur complement reduction to eliminate element unknowns and assemble the global system for face unknowns. Handles coupling terms for mixed systems and supports multi-dimensional problems.

  - uEquationHDG:
    High-level routine that calls the element, face, and Schur block assembly functions for all element blocks in the mesh.

  - RuEquationElemBlock:
    Computes only the residual vector for the element block, without assembling the full HDG matrices. Used for residual evaluation in nonlinear solvers.

  - RuEquationElemFaceBlock:
    Computes only the residual vector for the face block, including boundary and interface conditions.

  - ResidualHDG:
    High-level routine that computes the HDG residuals for all element blocks in the mesh.

  Notes:
  - The routines support both steady-state and time-dependent problems.
  - Auxiliary variables (e.g., w-equation) are supported for multiphysics coupling.
  - Boundary and interface conditions are imposed strongly via face contributions.
  - The code is designed for parallel execution and supports batched linear algebra operations.
  - Debugging output is available via file writing when debugMode is enabled.

  Dependencies:
  - Requires definitions for solstruct, resstruct, appstruct, masterstruct, meshstruct, tempstruct, commonstruct, and various linear algebra and utility routines.
  - Uses cublasHandle_t for GPU-accelerated batched matrix operations.

*/
#ifndef __UEQUATION
#define __UEQUATION

template <typename Model>
void uEquationElemBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int jth, Int backend)
{        
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int ncq = common.ncq;// number of compoments of (q)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg) 
    Int ncs = common.ncs;// number of compoments of (sdg) 
    Int ncw = common.ncw;// number of compoments of (wdg) 
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe; // number of nodes on master element
    Int nge = common.nge; // number of gauss points on master element        

    Int e1 = common.eblks[3*jth]-1;
    Int e2 = common.eblks[3*jth+1];            
    Int ne = e2-e1;
    //Int nn =  npe*ne; 
    Int nga = nge*ne;   
    //Int n0 = 0;                      // xg
    Int n1 = nga*ncx;                  // Xx
    Int n2 = nga*(ncx+nd*nd);          // jac        
    Int n3 = 0;                        // fg    
    Int n4 = nga*ncu*nd;               // sg  
    Int n5 = n4 + nga*ncu;             // ug
    Int n6 = n5 + nga*nc;              // wg
    Int n7 = n6 + nga*ncw;             // fg_uq
    Int n8 = n7 + nga*ncu*nd*nc;       // fg_w
    Int n9 = n8 + nga*ncu*nd*ncw;      // sg_uq
    Int n10 = n9 + nga*ncu*nc;         // sg_w
    Int n11 = n10 + nga*ncu*ncw;        // wg_uq
    Int nm = nge*e1*(ncx+nd*nd+1);

    // cout<<"j="<<jth<<endl;
    // cout<<"e1="<<e1<<endl;
    // cout<<"e2="<<e2<<endl;

    dstype *xg = &sol.elemg[nm];
    dstype *Xx = &sol.elemg[nm+n1];
    dstype *jac = &sol.elemg[nm+n2];

    dstype *og = &sol.odgg[nge*nco*e1];

    dstype *wsrc = &tmp.tempg[n3];        
    dstype *fg = &tmp.tempg[n3];        
    dstype *sg = &tmp.tempg[n4];
    dstype *uqg = &tmp.tempg[n5];
    dstype *wg = &tmp.tempg[n6];    
        
    dstype *fg_uq = &tmp.tempg[n7];    
    dstype *fg_w = &tmp.tempg[n8];    
    dstype *sg_uq = &tmp.tempg[n9];
    dstype *sg_w = &tmp.tempg[n10];    
    dstype *wg_uq = &tmp.tempg[n11];    

    // udg = tmp.tempg[n3] at gauss points on element
    GetElemNodes(tmp.tempn, sol.udg, npe, nc, 0, nc, e1, e2);
    //GetArrayAtIndex(tmp.tempn, sol.udg, &mesh.eindudg1[npe*nc*e1], nn*nc); // npe x ne x nc
    Node2Gauss(handle, uqg, tmp.tempn, master.shapegt, nge, npe, ne*nc, backend);    
    
    if ((ncw>0) & (common.wave==0)) {
        GetElemNodes(tmp.tempn, sol.wdg, npe, ncw, 0, ncw, e1, e2);    
        Node2Gauss(handle, wg, tmp.tempn, master.shapegt, nge, npe, ne*ncw, backend);        

        GetElemNodes(tmp.tempn, sol.wsrc, npe, ncw, 0, ncw, e1, e2);    
        Node2Gauss(handle, wsrc, tmp.tempn, master.shapegt, nge, npe, ne*ncw, backend);        
        
        // solve the w equation to get wg and wg_uq
        wEquation(wg, wg_uq, xg, uqg, og, wsrc, tmp.tempn, app, common, nga, backend);                
//         print2darray(uqg, nga, nc);
//         print2darray(wg, 1, nga);
//         print2darray(wg_uq, nga, nc);        
//         error("here");       
    }
    
    // calculate the source term Source(xdg, udg, odg, wdg)
    ArraySetValue(sg, 0.0, nga*ncu);
    ArraySetValue(sg_uq, 0.0, nga*ncu*nc);
    if ((ncw>0) & (common.wave==0)) ArraySetValue(sg_w, 0.0, nga*ncu*ncw);
    SourceDriver<Model>(sg, sg_uq, sg_w, xg, uqg, og, wg, mesh, master, app, sol, tmp, common, nge, e1, e2, backend);                 
    
    if (common.tdep) { // for time-dependent problem                
        // calculate sdg - udg*dtfactor
        ArrayAXPBY(fg, &sol.sdgg[nge*ncs*e1], uqg, one, -common.dtfactor, nga*ncu);                    
      
        if (common.tdfunc==1) {
            // calculate the time derivative function Tdfunc(xdg, udg, odg)
            TdfuncDriver<Model>(fg_uq, xg, uqg, og, wg, mesh, master, app, sol, tmp, common, nge, e1, e2, backend);            
        }
        else
            ArraySetValue(fg_uq, one, nga*ncu);;
                      
        // calculate (sdg-udg*dtfactor)*Tdfunc(xdg, udg, odg) 
        ArrayAXY(fg, fg, fg_uq, one, nga*ncu);                        
        
        // calculate sg = source + (sdgg - udg*dtfactor)*Tdfunc(xdg, udg, odg) 
        ArrayAXPBY(sg, sg, fg, one, one, nga*ncu);     
        
        // calculate sg_u = source_u  - dtfactor*Tdfunc(xdg, udg, odg) 
        ApplyDtcoef(sg_uq, fg_uq, -common.dtfactor, nga, ncu);                
        
//         // calculate sg = source + (sdgg - udg*dtfactor)*Tdfunc(xdg, udg, odg) 
//         ArrayAdd3Vectors(sg, sg, &sol.sdgg[nge*ncs*e1], uqg, one, one, -common.dtfactor, nga*ncu);        
//         
//         // calculate sg_u = source_u  - dtfactor*Tdfunc(xdg, udg, odg) 
//         ApplyDtcoef(sg_uq, -common.dtfactor, nga, ncu);                
    }
                
    // fg = nga*ncu*nd, sg = nga*ncu, fg_uq = nga*ncu*nd*nc, sg_uq = nga*ncu*nc
    ArraySetValue(fg, 0.0, nga*ncu*nd);
    ArraySetValue(fg_uq, 0.0, nga*ncu*nd*nc);
    if ((ncw>0) & (common.wave==0)) ArraySetValue(fg_w, 0.0, nga*ncu*nd*ncw); 
    FluxDriver<Model>(fg, fg_uq, fg_w, xg, uqg, og, wg, mesh, master, app, sol, tmp, common, nge, e1, e2, backend);    
    
    if ((ncw>0) & (common.wave==0)) {
        // sg_uq = sg_uq + sg_w * wg_uq -> ng * ncu * nc = ng * ncu * nc + (ng * ncu * ncw) * (ng * ncw * nc)
        ArrayGemmBatch2(sg_uq, sg_w, wg_uq, 1.0, ncu, nc, ncw, nga);
        // fg_uq = fg_uq + fg_w * wg_uq -> ng * ncu*nd * nc = ng * ncu*nd * nc + (ng * ncu*nd * ncw) * (ng * ncw * nc)
        ArrayGemmBatch2(fg_uq, fg_w, wg_uq, 1.0, ncu*nd, nc, ncw, nga);           
    }        
    
    // Ru = npe * ne * ncu 
    ApplyXxJac(tmp.tempn, sg, fg, Xx, jac, nge, nd, ncu, ne);
    Gauss2Node(handle, &res.Ru[npe*ncu*e1], tmp.tempn, master.shapegw, nge*(nd+1), npe, ncu*ne, backend);  // fixed bug here              

    // D = npe * npe * ne * ncu * ncu      
    ApplyXxJac(tmp.tempn, sg_uq, fg_uq, Xx, jac, nge, nd, ncu, ncu, ne);  // fixed bug here           
    Gauss2Node(handle, res.D, tmp.tempn, master.shapegwdotshapeg, nge*(nd+1), npe*npe, ncu*ncu*ne, backend);                
    ArrayMultiplyScalar(res.D, minusone, npe * npe * ne * ncu * ncu ); 

    if (ncq > 0) {
      // B = npe * npe * ne *  ncu * ncq 
      ApplyXxJac(tmp.tempn, &sg_uq[nga*ncu*ncu], &fg_uq[nga*ncu*nd*ncu], Xx, jac, nge, nd, ncu, ncq, ne);    // fixed bug here                    
      Gauss2Node(handle, res.B, tmp.tempn, master.shapegwdotshapeg, nge*(nd+1), npe*npe, ncu*ncq*ne, backend);        
      ArrayMultiplyScalar(res.B, minusone, npe * npe * ne * ncu * ncq );         
    }    

    if (common.debugMode==1) {      
      string filename;
      if (common.mpiProcs==1)
        filename = common.fileout;
      else
        filename = common.fileout + NumberToString(common.mpiRank);      

      writearray2file(filename + "uEquationElem_Ru.bin", res.Ru, npe*ncu*ne, backend);
      writearray2file(filename + "uEquationElem_D.bin", res.D, npe*npe*ncu*ncu*ne, backend);  
      if (ncq > 0) writearray2file(filename + "uEquationElem_B.bin", res.B, npe*npe*ncu*ncq*ne, backend);      
      writearray2file(filename + "uEquationElem_sg.bin", sg, nga*ncu, backend);
      writearray2file(filename + "uEquationElem_fg.bin", fg, nga*ncu*nd, backend);
      writearray2file(filename + "uEquationElem_sg_udg.bin", sg_uq, nga*ncu*nc, backend);
      writearray2file(filename + "uEquationElem_fg_udg.bin", fg_uq, nga*ncu*nd*nc, backend);
      writearray2file(filename + "uEquationElem_xg.bin", xg, nga*ncx, backend);
      writearray2file(filename + "uEquationElem_udg.bin", uqg, nga*nc, backend);
      if (ncw > 0) {
        writearray2file(filename + "uEquationElem_wg.bin", wg, nga*ncw, backend);
        writearray2file(filename + "uEquationElem_wg.bin", wg_uq, nga*ncw*nc, backend);
      }
    } 
}

template <typename Model>
void uEquationElemFaceBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int jth, Int backend)
{            
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int ncq = common.ncq;// number of compoments of (q)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int ncw = common.ncw;
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe; // number of nodes on master element
    Int npf = common.npf; // number of nodes on master face           
    Int ngf = common.ngf; // number of gauss poInts on master face              
    Int nfe = common.nfe; // number of faces in each element

    Int e1 = common.eblks[3*jth]-1;
    Int e2 = common.eblks[3*jth+1];            
    Int ne = e2-e1;
    Int nf = nfe*ne;
    Int nn =  npf*nf; 
    Int nga = ngf*nf;   
    Int n0 = 0;                                 // xg
    Int n1 = nga*ncx;                           // nlg
    Int n2 = nga*(ncx+nd);                      // jac
    //Int n3 = nga*(0);                           // uhg    
    Int n4 = nga*(ncu);                         // udg
    Int n5 = nga*(ncu+nc);                      // odg
    Int n6 = nga*(ncu+nc+nco);                  // wsrc
    Int n7 = nga*(ncu+nc+nco+ncw);              // wdg
    Int n8 = nga*(ncu+nc+nco+ncw+ncw);          // fhg
    Int nm = ngf*nfe*e1*(ncx+nd+1);
    
    dstype *xg = &sol.elemfaceg[nm+n0];    
    dstype *nlg = &sol.elemfaceg[nm+n1];
    dstype *jac = &sol.elemfaceg[nm+n2];

    dstype *uhg = &tmp.tempg[0];
    dstype *udg = &tmp.tempg[n4];
    dstype *odg = &tmp.tempg[n5];
    dstype *wsrc = &tmp.tempg[n6];
    dstype *wdg = &tmp.tempg[n7];

    dstype *fh     = &tmp.tempg[n8];
    dstype *fh_uq  = &tmp.tempg[n8 + nga*ncu*nd];
    dstype *fh_uh  = &tmp.tempg[n8 + nga*ncu*nd + nga*ncu*nd*nc];
    dstype *fh_w   = &tmp.tempg[n8 + nga*ncu*nd + nga*ncu*nd*nc + nga*ncu*ncu];
    dstype *wdg_uq = &tmp.tempg[n8 + nga*ncu*nd + nga*ncu*nd*nc + nga*ncu*ncu + nga*ncu*nd*ncw];        
                
    // npf * nfe * ne * ncu
    GetElementFaceNodes(tmp.tempn, sol.uh, mesh.elemcon, npf*nfe, ncu, e1, e2, 0); // fixed bug here

    // udg = tmp.tempg[n4] at gauss points on face
    GetElementFaceNodes(&tmp.tempn[nn*ncu], sol.udg, mesh.perm, npf*nfe, nc, npe, nc, e1, e2);

    if (nco>0) GetElementFaceNodes(&tmp.tempn[nn*(ncu+nc)], sol.odg, mesh.perm, npf*nfe, nco, npe, nco, e1, e2);      

    if ((ncw>0) & (common.wave==0)) {
      GetElementFaceNodes(&tmp.tempn[nn*(ncu+nc+nco)], sol.wsrc, mesh.perm, npf*nfe, ncw, npe, ncw, e1, e2); 
      GetElementFaceNodes(&tmp.tempn[nn*(ncu+nc+nco+ncw)], sol.wdg, mesh.perm, npf*nfe, ncw, npe, ncw, e1, e2); // fix bug here
    }
    
    Node2Gauss(handle, tmp.tempg, tmp.tempn, master.shapfgt, ngf, npf, nfe*ne*(ncu+nc+nco+ncw+ncw), backend); // fix bug here
        
    if ((ncw>0) & (common.wave==0)) {
        // copy udg to tmp.tempn
        ArrayCopy(tmp.tempn, udg, nga*nc);
        
        // replace u with uhat 
        ArrayCopy(tmp.tempn, uhg, nga*ncu);
            
        // solve the w equation to get wg and wg_uq
        wEquation(wdg, wdg_uq, xg, tmp.tempn, odg, wsrc, &tmp.tempn[nga*nc], app, common, nga, backend);
                
//         print2darray(wdg, ngf*nfe, ncw, nga, ncw);
//         error("here");
        
        // solve the w equation to get wg and wg_uq
        // wEquation(wdg, wdg_uq, xg, udg, odg, wsrc, tmp.tempn, app, common, nga, backend);
    }

    ArraySetValue(fh, 0.0, nga*ncu);    
    ArraySetValue(fh_uq, 0.0, nga*ncu*nc);        
    ArraySetValue(fh_uh, 0.0, nga*ncu*ncu);
    
    if (ncw > 0) ArraySetValue(fh_w, 0.0, nga*ncu*ncw);       
            
    FhatDriver<Model>(fh, fh_uq, fh_w, fh_uh, xg, udg, odg, wdg, uhg, nlg, 
        mesh, master, app, sol, tmp, common, nga, backend);      
        
    if ((ncw>0) & (common.wave==0)) {
      ArrayGemmBatch2(fh_uh, fh_w, wdg_uq, one, ncu, ncu, ncw, nga); // fix bug here       
      
      ArraySetValue(wdg_uq, 0.0, nga*ncu*ncu);
      ArrayGemmBatch2(fh_uq, fh_w, wdg_uq, one, ncu, nc, ncw, nga); // fix bug here       
    }
    
//     if (common.debugMode==1) {    
//       writearray2file(common.fileout + "uEquationElemFace_fh.bin", fh, nga*ncu, backend);
//       writearray2file(common.fileout + "uEquationElemFace_fh_udg.bin", fh_uq, nga*ncu, backend);
//     }
    
    //columnwiseMultiply(dstype* C, const dstype* A, const dstype* b, const int N, const int M)
    columnwiseMultiply(fh, fh, jac, nga, ncu);
    columnwiseMultiply(fh_uq, fh_uq, jac, nga, ncu*nc);
    columnwiseMultiply(fh_uh, fh_uh, jac, nga, ncu*ncu);

    dstype *Rutmp = &tmp.tempn[0];
    dstype *Dtmp  = &tmp.tempn[npf*nfe*ne*ncu];
    dstype *Btmp  = &tmp.tempn[npf*nfe*ne*ncu + npf*npf*nfe*ne*ncu*ncu];
    dstype *Ftmp  = &tmp.tempn[npf*nfe*ne*ncu + npf*npf*nfe*ne*ncu*ncu + npf*npf*nfe*ne*ncu*ncq];
    
    // npf*nfe*ne*ncu
    Gauss2Node(handle, Rutmp, fh, master.shapfgw, ngf, npf, nf*ncu, backend);                

    // if (common.debugMode==1) {    
    //   writearray2file(common.fileout + "uEquationElemFace_Rutmp.bin", Rutmp, npf*nfe*ne*ncu, backend);
    // }

    // npf*nfe*ne*ncu -> npe*ne*ncu
    assembleRu(&res.Ru[npe*ncu*e1], Rutmp, mesh.perm, npe, npf*nfe, ne*ncu);

    // npf*npf*nfe*ne*ncu*ncu
    Gauss2Node(handle, Dtmp, fh_uq, master.shapfgwdotshapfg, ngf, npf*npf, nf*ncu*ncu, backend);            
    
    // npf*npf*nfe*ne*ncu*ncu -> npe*npe*ne*ncu*ncu  
    assembleMatrixBD(res.D, Dtmp, mesh.perm, npe, npf, nfe, ne*ncu*ncu);

    if (ncq > 0) {
      // npf*npf*nfe*ne*ncu*ncq
      Gauss2Node(handle, Btmp, &fh_uq[nga*ncu*ncu], master.shapfgwdotshapfg, ngf, npf*npf, nf*ncu*ncq, backend);            
      // npf*npf*nfe*ne*ncu*ncq -> npe*npe*ne*ncu*ncq
      assembleMatrixBD(res.B, Btmp, mesh.perm, npe, npf, nfe, ne*ncu*ncq);
    }

    // npf*npf*nfe*ne*ncu*ncu
    Gauss2Node(handle, Ftmp, fh_uh, master.shapfgwdotshapfg, ngf, npf*npf, nf*ncu*ncu, backend);            
    // npf*npf*nfe*ne*ncu*ncu -> npe*npf*nfe*ne*ncu*ncu
    assembleMatrixF(&res.F[npe*npf*nfe*ncu*ncu*e1], Ftmp, mesh.perm, npe, npf, nfe, ne*ncu*ncu);

    // impose boundary conditions
    for (int ibc=0; ibc<common.maxnbc; ibc++)
    {
      int n = ibc + common.maxnbc*jth;
      int start = common.nboufaces[n];
      int nfaces = common.nboufaces[n + 1] - start;
      if (nfaces>0) {        
        int ngb = nfaces*ngf;
        dstype *xgb = &tmp.tempg[n8];
        dstype *ugb = &tmp.tempg[n8 + ngb*ncx];
        dstype *ogb = &tmp.tempg[n8 + ngb*ncx + ngb*nc];
        dstype *wgb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco]; 
        dstype *uhb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw];
        dstype *nlb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu];
        dstype *wsb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu + ngb*nd];
        dstype *fhb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu + ngb*nd + ngb*ncw];
        dstype *fhb_uq = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu + ngb*nd + ngb*ncw + ngb*ncu];
        dstype *fhb_w = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu + ngb*nd + ngb*ncw + ngb*ncu + ngb*ncu*nc];
        dstype *fhb_uh = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu + ngb*nd + ngb*ncw + ngb*ncu + ngb*ncu*nc + ngb*ncu*ncw];
        dstype *wgb_uq = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu + ngb*nd + ngb*ncw + ngb*ncu + ngb*ncu*nc + ngb*ncu*ncw + ngb*ncu*ncu];
        //dstype *Rb =  &tmp.tempn[npf*nfe*ne*ncu + npf*npf*nfe*ne*ncu*ncu + npf*npf*nfe*ne*ncu*ncq + npf*npf*nfe*ne*ncu*ncu];
                
        GetBoundaryNodes(xgb, xg, &mesh.boufaces[start], ngf, nfe, ne, ncx, nfaces);
        GetBoundaryNodes(ugb, udg, &mesh.boufaces[start], ngf, nfe, ne, nc, nfaces);
        GetBoundaryNodes(ogb, odg, &mesh.boufaces[start], ngf, nfe, ne, nco, nfaces);
        GetBoundaryNodes(wgb, wdg, &mesh.boufaces[start], ngf, nfe, ne, ncw, nfaces);
        GetBoundaryNodes(wsb, wsrc, &mesh.boufaces[start], ngf, nfe, ne, ncw, nfaces);
        GetBoundaryNodes(uhb, uhg, &mesh.boufaces[start], ngf, nfe, ne, ncu, nfaces);
        GetBoundaryNodes(nlb, nlg, &mesh.boufaces[start], ngf, nfe, ne, nd, nfaces);

        if ((ncw>0) & (common.wave==0)) {
          // copy (u, q) to res.K
          ArrayCopy(res.K, ugb, ngb*nc);
        
          // replace u with uhat 
          ArrayCopy(res.K, uhb, ngb*ncu);
        
          wEquation(wgb, wgb_uq, xgb, res.K, ogb, wsb, &res.K[ngb*nc], app, common, ngb, backend);          
          
          // wEquation(wgb, wgb_uq, xgb, ugb, ogb, wsb, Rb, app, common, ngb, backend);
        }
        
        // intialize fhb, fhb_uq, fhb_w, fhb_uh to zero 
        ArraySetValue(fhb, 0.0, ngb*ncu);
        ArraySetValue(fhb_uq, 0.0, ngb*ncu*nc);
        ArraySetValue(fhb_uh, 0.0, ngb*ncu*ncu);
        if (ncw > 0) ArraySetValue(fhb_w, 0.0, ngb*ncu*ncw);      
        if (ibc+1 == 1000) 
            StgInflowHDG(fhb, fhb_uq, fhb_w, fhb_uh, res.K, xgb, ogb, uhb, 
                         app.physicsparam, app.stgdata, app.stgparam, common.time, 
                         ngb, common.stgNmode, nd, ncu, nc, ncw);
         else
            FbouDriver<Model>(fhb, fhb_uq, fhb_w, fhb_uh, xgb, ugb, ogb, wgb, uhb, nlb, 
                 mesh, master, app, sol, tmp, common, ngb, ibc+1, backend);    

        if ((ncw>0) & (common.wave==0)) {      
          // void ArrayGemmBatch2(dstype* C, const dstype* A, const dstype* B, dstype alpha, const int I, const int J, const int K, const int S)
          // C[S*I*J] = A[S*I*K] x B[S*K*J] + C[S*I*J]
          ArrayGemmBatch2(fhb_uh, fhb_w, wgb_uq, one, ncu, ncu, ncw, ngb);  // fix bug here             
          //ArraySetValue(wgb_uq, 0.0, ngb*ncu*ncu);          
          ArraySetValue(wgb_uq, 0.0, ngb*ncw*ncu);          
          ArrayGemmBatch2(fhb_uq, fhb_w, wgb_uq, one, ncu, nc, ncw, ngb);  // fix bug here            
          //ArrayGemmBatch2(fhb_uq, fhb_w, wgb_uq, one, ncu, nc, ncw, ngb);  // fix bug here                 
        }

        dstype *jacb = &tmp.tempg[n8];
        GetBoundaryNodes(jacb, jac, &mesh.boufaces[start], ngf, nfe, ne, 1, nfaces);
        columnwiseMultiply(fhb, fhb, jacb, ngb, ncu);
        columnwiseMultiply(fhb_uq, fhb_uq, jacb, ngb, ncu*nc);
        columnwiseMultiply(fhb_uh, fhb_uh, jacb, ngb, ncu*ncu);

        dstype *Rb =  res.K;
        Gauss2Node(handle, Rb, fhb, master.shapfgw, ngf, npf, nfaces*ncu, backend);                
        PutBoundaryNodes(Rutmp, Rb, &mesh.boufaces[start], npf, nfe, ne, ncu, nfaces);        

        Gauss2Node(handle, Rb, fhb_uq, master.shapfgwdotshapfg, ngf, npf*npf, nfaces*ncu*ncu, backend);
        PutBoundaryNodes(Dtmp, Rb, &mesh.boufaces[start], npf*npf, nfe, ne, ncu*ncu, nfaces);

        if (ncq > 0) {
          Gauss2Node(handle, Rb, &fhb_uq[ngb*ncu*ncu], master.shapfgwdotshapfg, ngf, npf*npf, nfaces*ncu*ncq, backend);          
          PutBoundaryNodes(Btmp, Rb, &mesh.boufaces[start], npf*npf, nfe, ne, ncu*ncq, nfaces);
        }

        Gauss2Node(handle, Rb, fhb_uh, master.shapfgwdotshapfg, ngf, npf*npf, nfaces*ncu*ncu, backend);
        PutBoundaryNodes(Ftmp, Rb, &mesh.boufaces[start], npf*npf, nfe, ne, ncu*ncu, nfaces);                
      }
    }

    // npf*nfe*ne*ncu -> npf*nfe*ne*ncu
    ArrayAXPB(&res.Rh[npf*nfe*ncu*e1], Rutmp, minusone, zero, npf*nfe*ne*ncu);

    // npf*npf*nfe*ne*ncu*ncu -> npf*nfe*npe*ne*ncu*ncu
    ArraySetValue(res.K, zero, npf*nfe*npe*ne*ncu*ncu);
    assembleMatrixGK(res.K, Dtmp, mesh.perm, npe, npf, nfe, ne*ncu*ncu);

    if (ncq > 0) {
      // npf*npf*nfe*ne*ncu*ncq -> npf*nfe*npe*ne*ncu*ncq
      ArraySetValue(res.G, zero, npf*nfe*npe*ne*ncu*ncq);
      assembleMatrixGK(res.G, Btmp, mesh.perm, npe, npf, nfe, ne*ncu*ncq);
    }

    // npf*npf*nfe*ne*ncu*ncu -> npf*nfe*npf*nfe*ne*ncu*ncu
    assembleMatrixH(&res.H[npf*nfe*npf*nfe*ncu*ncu*e1], Ftmp, mesh.perm, npe, npf, nfe, ne*ncu*ncu);

    // impose interface conditions
    for (int ibc=0; ibc<common.maxnbc; ibc++)
    {
      int n = ibc + common.maxnbc*jth;
      int start = common.nboufaces[n];
      int nfaces = common.nboufaces[n + 1] - start;
      if ((nfaces>0) && (common.eblks[3*jth+2]==-1) && (common.coupledboundarycondition == ibc+1) && (common.coupledcondition>0)) {        
        int ncu12 = common.szinterfacefluxmap;
        
        //printf("%d %d %d %d %d %d %d %d %d %d %d\n", common.mpiRank, common.eblks[3*jth], common.eblks[3*jth+1], common.eblks[3*jth+2], ibc, common.coupledinterface, common.coupledcondition, common.coupledboundarycondition, common.ncie, ne, ncu12);
                        
        int ngb = nfaces*ngf;
        dstype *xgb = &tmp.tempg[n8];
        dstype *ugb = &tmp.tempg[n8 + ngb*ncx];
        dstype *ogb = &tmp.tempg[n8 + ngb*ncx + ngb*nc];
        dstype *wgb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco]; 
        dstype *uhb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw];
        dstype *nlb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu];
        dstype *wsb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu + ngb*nd];
        dstype *fhb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu + ngb*nd + ngb*ncw];
        dstype *fhb_uq = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu + ngb*nd + ngb*ncw + ngb*ncu];
        dstype *fhb_w = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu + ngb*nd + ngb*ncw + ngb*ncu + ngb*ncu*nc];
        dstype *fhb_uh = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu + ngb*nd + ngb*ncw + ngb*ncu + ngb*ncu*nc + ngb*ncu*ncw];
        dstype *wgb_uq = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu + ngb*nd + ngb*ncw + ngb*ncu + ngb*ncu*nc + ngb*ncu*ncw + ngb*ncu*ncu];               
        
        GetBoundaryNodes(xgb, xg, &mesh.boufaces[start], ngf, nfe, ne, ncx, nfaces);
        GetBoundaryNodes(ugb, udg, &mesh.boufaces[start], ngf, nfe, ne, nc, nfaces);
        GetBoundaryNodes(ogb, odg, &mesh.boufaces[start], ngf, nfe, ne, nco, nfaces);
        GetBoundaryNodes(wgb, wdg, &mesh.boufaces[start], ngf, nfe, ne, ncw, nfaces);
        GetBoundaryNodes(wsb, wsrc, &mesh.boufaces[start], ngf, nfe, ne, ncw, nfaces);
        GetBoundaryNodes(uhb, uhg, &mesh.boufaces[start], ngf, nfe, ne, ncu, nfaces);
        GetBoundaryNodes(nlb, nlg, &mesh.boufaces[start], ngf, nfe, ne, nd, nfaces);

        // transfer summit's boundary data to ogb  

        if ((ncw>0) & (common.wave==0)) {
          dstype *temp1 =  &tmp.tempn[0];
          dstype *temp2 =  &tmp.tempn[ngb*nc];
          
          // copy (u, q) to temp1
          ArrayCopy(temp1, ugb, ngb*nc);
        
          // replace u with uhat 
          ArrayCopy(temp1, uhb, ngb*ncu);
        
          wEquation(wgb, wgb_uq, xgb, temp1, ogb, wsb, temp2, app, common, ngb, backend);                    
        }                
                
        // intialize fhb, fhb_uq, fhb_w, fhb_uh to zero 
        ArraySetValue(fhb, 0.0, ngb*ncu12);
        ArraySetValue(fhb_uq, 0.0, ngb*ncu12*nc);
        ArraySetValue(fhb_uh, 0.0, ngb*ncu12*ncu);
        if (ncw > 0) ArraySetValue(fhb_w, 0.0, ngb*ncu12*ncw);        
        FintDriver<Model>(fhb, fhb_uq, fhb_w, fhb_uh, xgb, ugb, ogb, wgb, uhb, nlb, 
             mesh, master, app, sol, tmp, common, ngb, common.coupledcondition, backend);    
                        
        if ((ncw>0) & (common.wave==0)) {          
          ArrayGemmBatch2(fhb_uh, fhb_w, wgb_uq, one, ncu12, ncu, ncw, ngb);  // fix bug here             
          ArraySetValue(wgb_uq, 0.0, ngb*ncw*ncu);          
          ArrayGemmBatch2(fhb_uq, fhb_w, wgb_uq, one, ncu12, nc, ncw, ngb);  // fix bug here                      
        }

        dstype *jacb = &tmp.tempg[n8];
        GetBoundaryNodes(jacb, jac, &mesh.boufaces[start], ngf, nfe, ne, 1, nfaces);
        columnwiseMultiply(fhb, fhb, jacb, ngb, ncu12);
        columnwiseMultiply(fhb_uq, fhb_uq, jacb, ngb, ncu12*nc);
        columnwiseMultiply(fhb_uh, fhb_uh, jacb, ngb, ncu12*ncu);

        dstype *Rb =  &tmp.tempn[0];
        
        Gauss2Node(handle, Rb, fhb, master.shapfgw, ngf, npf, nfaces*ncu12, backend);                
        ArrayAXPB(res.Ri, Rb, minusone, zero, npf*nfaces*ncu12);
        
        ArraySetValue(res.Ki, zero, ncu12*npf*npe*ncu*nfaces);
        Gauss2Node(handle, Rb, fhb_uq, master.shapfgwdotshapfg, ngf, npf*npf, nfaces*ncu12*ncu, backend);        
        // npf*npf*nfaces*ncu12*ncu -> ncu12*npf*npe*ncu*nfaces        
        assembleMatrixKint(res.Ki, Rb, &mesh.boufaces[start], mesh.perm, npe, npf, nfe, ncu12, ncu, nfaces);        
                
        
        if (ncq > 0) {
          ArraySetValue(res.Gi, zero, npf*npe*nfaces*ncu12*ncq);
          Gauss2Node(handle, Rb, &fhb_uq[ngb*ncu12*ncu], master.shapfgwdotshapfg, ngf, npf*npf, nfaces*ncu12*ncq, backend);          
          // npf*npf*nfaces*ncu12*ncq ->  npf*npe*nfaces*ncu12*ncq
          assembleMatrixGint(res.Gi, Rb, &mesh.boufaces[start], mesh.perm, npe, npf, nfe, ncu12, ncq, nfaces);
        }
        
        ArraySetValue(res.Hi, zero, ncu12*npf*ncu*npf*nfe*nfaces);
        Gauss2Node(handle, Rb, fhb_uh, master.shapfgwdotshapfg, ngf, npf*npf, nfaces*ncu12*ncu, backend);
        // npf*npf*nfaces*ncu12*ncu -> ncu12*npf*ncu*npf*nfe*nfaces        
        assembleMatrixHint(res.Hi, Rb, &mesh.boufaces[start], npe, npf, nfe, ncu12, ncu, nfaces);                  
      }      
    }
    
  if (common.debugMode==1) {    
    string filename;
    if (common.mpiProcs==1)
      filename = common.fileout;
    else
      filename = common.fileout + NumberToString(common.mpiRank);      
    writearray2file(filename + "uEquationElemFace_Ru.bin", &res.Ru[npe*ncu*e1], npe*ne*ncu, backend);
    writearray2file(filename + "uEquationElemFace_D.bin", res.D, npe*npe*ncu*ncu*ne, backend);
    writearray2file(filename + "uEquationElemFace_B.bin", res.B, npe*npe*ncu*ncq*ne, backend);
    writearray2file(filename + "uEquationElemFace_F.bin", &res.F[npe*npf*nfe*ncu*ncu*e1], npe*npf*nfe*ncu*ncu*ne, backend);
    writearray2file(filename + "uEquationElemFace_Rh.bin", &res.Rh[npf*nfe*ncu*e1], npf*nfe*ne*ncu, backend);
    writearray2file(filename + "uEquationElemFace_K.bin", res.K, npf*nfe*npe*ne*ncu*ncu, backend);  
    writearray2file(filename + "uEquationElemFace_G.bin", res.G, npf*nfe*npe*ne*ncu*ncq, backend);  
    writearray2file(filename + "uEquationElemFace_H.bin", &res.H[npf*nfe*npf*nfe*ncu*ncu*e1], npf*nfe*npf*nfe*ne*ncu*ncu, backend);      
  }
}

void uEquationSchurBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int jth, Int backend)
{        
    Int ncu = common.ncu;// number of compoments of (u)
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe; // number of nodes on master element
    Int npf = common.npf; // number of nodes on master face           
    Int nfe = common.nfe; // number of faces in each element

    Int e1 = common.eblks[3*jth]-1;
    Int e2 = common.eblks[3*jth+1];            
    Int ne = e2-e1;

    Int n = npe*ncu; 
    Int m = npf*nfe*ncu;
    dstype *DinvF = &res.F[n*m*e1];    
    dstype *Ru = &res.Ru[n*e1];
    dstype *DinvH = &res.H[m*m*e1];
    dstype *Rh = &res.Rh[m*e1];

//     dstype *D = tmp.tempn;
//     dstype *F = &tmp.tempn[npe*ncu*npe*ncu*ne];
//     dstype *K = &tmp.tempn[npe*ncu*npe*ncu*ne + npe*npf*nfe*ncu*ncu*ne];
//     dstype *H = &tmp.tempn[npe*ncu*npe*ncu*ne + npe*npf*nfe*ncu*ncu*ne + npe*npf*nfe*ncu*ncu*ne];
// 
//     // npe * npe * ne * ncu * ncu -> npe * ncu * npe * ncu * ne
//     schurMatrixD(D, res.D, npe, ncu, ne);
//     // npe * npf * nfe * ne * ncu * ncu -> npe * ncu * ncu * npf * nfe * ne
//     schurMatrixF(F, &res.F[n*m*e1], npe, ncu, npf, nfe, ne);
//     // npf * nfe * npe * ne * ncu * ncu -> ncu * npf * nfe * npe * ncu * ne
//     schurMatrixK(K, res.K, npe, ncu, ncu, npf, nfe, ne);    
//     // npf * nfe * npf * nfe * ne * ncu * ncu -> ncu * npf * nfe * ncu * npf * nfe * ne
//     schurMatrixH(H, &res.H[m*m*e1], ncu, ncu, npf, nfe, ne);

    dstype *D = res.D;
    dstype *K = res.K;
    dstype *F = &res.F[n*m*e1];
    dstype *H = &res.H[m*m*e1];

    // npe * npe * ne * ncu * ncu -> npe * ncu * npe * ncu * ne
    schurMatrixD(tmp.tempn, res.D, npe, ncu, ne);
    ArrayCopy(D, tmp.tempn, n * n * ne);
            
    // npf * nfe * npe * ne * ncu * ncu -> ncu * npf * nfe * npe * ncu * ne
    schurMatrixK(tmp.tempn, res.K, npe, ncu, ncu, npf, nfe, ne);    
    ArrayCopy(K, tmp.tempn, n * m * ne);
    
    // npf * nfe * npf * nfe * ne * ncu * ncu -> ncu * npf * nfe * ncu * npf * nfe * ne
    schurMatrixH(tmp.tempn, &res.H[m*m*e1], ncu, ncu, npf, nfe, ne);
    ArrayCopy(H, tmp.tempn, m * m * ne);
    
    // npe * npf * nfe * ne * ncu * ncu -> npe * ncu * ncu * npf * nfe * ne
    schurMatrixF(tmp.tempn, &res.F[n*m*e1], npe, ncu, npf, nfe, ne);
    ArrayCopy(F, tmp.tempn, n * m * ne);
        
    dstype scalar = 1.0;
    if (common.wave==1)
        scalar = 1.0/common.dtfactor;    
    
    if (common.ncq > 0) {      
      if (nd == 1) {
        // D = D + B * Minv * C
        schurMatrixBMinvC(D, res.B, &res.C[npe*npe*e1], scalar, npe, ncu, ne);
        // F = F - B * Minv * E
        schurMatrixBMinvE(F, res.B, &res.E[npe*npf*nfe*e1], scalar, npe, ncu, npf, nfe, ne);
        // K = K + G * Minv * C
        schurMatrixGMinvC(K, res.G, &res.C[npe*npe*e1], scalar, npe, ncu, ncu, npf, nfe, ne);
        // H = H - G * Minv * E
        schurMatrixGMinvE(H, res.G, &res.E[npe*npf*nfe*e1], scalar, npe, ncu, ncu, npf, nfe, ne);
      } 
      else if (nd == 2) {
        dstype *Cx = &res.C[npe*npe*e1]; // fix bug here
        dstype *Cy = &res.C[npe*npe*common.ne + npe*npe*e1]; // fix bug here
        dstype *Ex = &res.E[npe*npf*nfe*e1]; // fix bug here
        dstype *Ey = &res.E[npe*npf*nfe*common.ne + npe*npf*nfe*e1]; // fix bug here
        dstype *Bx = res.B; // npe * npe * ne * ncu * ncu
        dstype *By = &res.B[npe*npe*ncu*ncu*ne]; // npe * npe * ne * ncu * ncu
        dstype *Gx = res.G; // npf*nfe*npe*ne*ncu*ncu
        dstype *Gy = &res.G[npf*nfe*npe*ncu*ncu*ne]; // npf*nfe*npe*ne*ncu*ncu

        schurMatrixBMinvC(D, Bx, Cx, scalar, npe, ncu, ne);
        schurMatrixBMinvC(D, By, Cy, scalar, npe, ncu, ne);

        schurMatrixBMinvE(F, Bx, Ex, scalar, npe, ncu, npf, nfe, ne);
        schurMatrixBMinvE(F, By, Ey, scalar, npe, ncu, npf, nfe, ne);

        schurMatrixGMinvC(K, Gx, Cx, scalar, npe, ncu, ncu, npf, nfe, ne);
        schurMatrixGMinvC(K, Gy, Cy, scalar, npe, ncu, ncu, npf, nfe, ne);

        schurMatrixGMinvE(H, Gx, Ex, scalar, npe, ncu, ncu, npf, nfe, ne); 
        schurMatrixGMinvE(H, Gy, Ey, scalar, npe, ncu, ncu, npf, nfe, ne);
      }
      else if (nd == 3) {
        dstype *Cx = &res.C[npe*npe*e1]; // fixed bug here
        dstype *Cy = &res.C[npe*npe*common.ne + npe*npe*e1]; // fixed bug here
        dstype *Cz = &res.C[npe*npe*common.ne*2 + npe*npe*e1]; // fixed bug here
        dstype *Ex = &res.E[npe*npf*nfe*e1]; // fixed bug here
        dstype *Ey = &res.E[npe*npf*nfe*common.ne + npe*npf*nfe*e1]; // fixed bug here
        dstype *Ez = &res.E[npe*npf*nfe*common.ne*2 + npe*npf*nfe*e1]; // fixed bug here
        dstype *Bx = res.B;
        dstype *By = &res.B[npe*npe*ncu*ncu*ne];
        dstype *Bz = &res.B[npe*npe*ncu*ncu*ne*2];
        dstype *Gx = res.G;
        dstype *Gy = &res.G[npf*nfe*npe*ncu*ncu*ne];
        dstype *Gz = &res.G[npf*nfe*npe*ncu*ncu*ne*2];

        schurMatrixBMinvC(D, Bx, Cx, scalar, npe, ncu, ne);
        schurMatrixBMinvC(D, By, Cy, scalar, npe, ncu, ne);
        schurMatrixBMinvC(D, Bz, Cz, scalar, npe, ncu, ne);

        schurMatrixBMinvE(F, Bx, Ex, scalar, npe, ncu, npf, nfe, ne);
        schurMatrixBMinvE(F, By, Ey, scalar, npe, ncu, npf, nfe, ne);
        schurMatrixBMinvE(F, Bz, Ez, scalar, npe, ncu, npf, nfe, ne);

        schurMatrixGMinvC(K, Gx, Cx, scalar, npe, ncu, ncu, npf, nfe, ne);
        schurMatrixGMinvC(K, Gy, Cy, scalar, npe, ncu, ncu, npf, nfe, ne);
        schurMatrixGMinvC(K, Gz, Cz, scalar, npe, ncu, ncu, npf, nfe, ne);

        schurMatrixGMinvE(H, Gx, Ex, scalar, npe, ncu, ncu, npf, nfe, ne);
        schurMatrixGMinvE(H, Gy, Ey, scalar, npe, ncu, ncu, npf, nfe, ne);
        schurMatrixGMinvE(H, Gz, Ez, scalar, npe, ncu, ncu, npf, nfe, ne);
      }
    }

    if (common.debugMode==1) {    
      string filename;
      if (common.mpiProcs==1)
        filename = common.fileout;
      else
        filename = common.fileout + NumberToString(common.mpiRank);              
      writearray2file(filename + "uEquationElemSchur_D.bin", D, npe*npe*ncu*ncu*ne, backend);
      writearray2file(filename + "uEquationElemSchur_F.bin", F, npe*npf*nfe*ncu*ncu*ne, backend);
      writearray2file(filename + "uEquationElemSchur_K.bin", K, npe*npf*nfe*ne*ncu*ncu, backend);  
      writearray2file(filename + "uEquationElemSchur_H.bin", H, npf*npf*nfe*nfe*ne*ncu*ncu, backend);      
    } 

    // compute the inverse of D 
    Inverse(handle, D, tmp.tempn, res.ipiv, n, ne, backend);    
    
    // DinvF = -Dinv * F
    ArrayCopy(tmp.tempn, F, n * m * ne);
    PGEMNMStridedBached(handle, n, m, n, minusone, D, n, tmp.tempn, n, zero, DinvF, n, ne, backend); // fixed bug here

    // Ru = Dinv * Ru    
    schurVectorRu(tmp.tempn, Ru, npe, ncu, ne); // permute Ru from npe*ne*ncu to npe*ncu*ne
    PGEMNMStridedBached(handle, n, 1, n, one, D, n, tmp.tempn, n, zero, Ru, n, ne, backend);    

    // H = H + K * DinvF
    //ArrayCopy(DinvH, H, m*m*ne);
    PGEMNMStridedBached(handle, m, m, n, one, K, m, DinvF, n, one, DinvH, m, ne, backend);

    // Rh = Rh - K * Dinv * Ru
    schurVectorRh(tmp.tempn, Rh, npf*nfe, ncu, ne); // permute Rh from npf*nfe*ne*ncu to ncu*npf*nfe*ne
    ArrayCopy(Rh, tmp.tempn, m*ne);
    PGEMNMStridedBached(handle, m, 1, n, minusone, K, m, Ru, n, one, Rh, m, ne, backend);     
    
    if ((common.eblks[3*jth+2]==-1) && (common.coupledcondition>0) && (common.coupledinterface>0)) {      
      int ncu12 = common.szinterfacefluxmap;
      Int m12 = npf*ncu12;
      
      if (common.ncq > 0) {      
        if (nd == 1) {
          schurMatrixGMinvC(res.Ki, res.Gi, &res.C[npe*npe*e1], scalar, npe, ncu12, ncu, npf, 1, ne);
          schurMatrixGintMinvE(res.Hi, res.Gi, &res.E[npe*npf*nfe*e1], scalar, npe, ncu12, ncu, npf, nfe, ne);
        } 
        else if (nd == 2) {
          dstype *Cx = &res.C[npe*npe*e1]; // fix bug here
          dstype *Cy = &res.C[npe*npe*common.ne + npe*npe*e1]; // fix bug here
          dstype *Ex = &res.E[npe*npf*nfe*e1]; // fix bug here
          dstype *Ey = &res.E[npe*npf*nfe*common.ne + npe*npf*nfe*e1]; // fix bug here
          dstype *Gx = res.Gi; // npf*npe*ne*ncu12*ncu
          dstype *Gy = &res.Gi[npf*npe*ncu12*ncu*ne]; // npf*npe*ne*ncu12*ncu
          
          schurMatrixGMinvC(res.Ki, Gx, Cx, scalar, npe, ncu12, ncu, npf, 1, ne);
          schurMatrixGMinvC(res.Ki, Gy, Cy, scalar, npe, ncu12, ncu, npf, 1, ne);
          schurMatrixGintMinvE(res.Hi, Gx, Ex, scalar, npe, ncu12, ncu, npf, nfe, ne); 
          schurMatrixGintMinvE(res.Hi, Gy, Ey, scalar, npe, ncu12, ncu, npf, nfe, ne);
        }
        else if (nd == 3) {
          dstype *Cx = &res.C[npe*npe*e1]; // fixed bug here
          dstype *Cy = &res.C[npe*npe*common.ne + npe*npe*e1]; // fixed bug here
          dstype *Cz = &res.C[npe*npe*common.ne*2 + npe*npe*e1]; // fixed bug here
          dstype *Ex = &res.E[npe*npf*nfe*e1]; // fixed bug here
          dstype *Ey = &res.E[npe*npf*nfe*common.ne + npe*npf*nfe*e1]; // fixed bug here
          dstype *Ez = &res.E[npe*npf*nfe*common.ne*2 + npe*npf*nfe*e1]; // fixed bug here
          dstype *Gx = res.Gi;
          dstype *Gy = &res.Gi[npf*npe*ncu12*ncu*ne];
          dstype *Gz = &res.Gi[npf*npe*ncu12*ncu*ne*2];

          schurMatrixGMinvC(res.Ki, Gx, Cx, scalar, npe, ncu12, ncu, npf, 1, ne);
          schurMatrixGMinvC(res.Ki, Gy, Cy, scalar, npe, ncu12, ncu, npf, 1, ne);
          schurMatrixGMinvC(res.Ki, Gz, Cz, scalar, npe, ncu12, ncu, npf, 1, ne);
          schurMatrixGintMinvE(res.Hi, Gx, Ex, scalar, npe, ncu12, ncu, npf, nfe, ne);
          schurMatrixGintMinvE(res.Hi, Gy, Ey, scalar, npe, ncu12, ncu, npf, nfe, ne);
          schurMatrixGintMinvE(res.Hi, Gz, Ez, scalar, npe, ncu12, ncu, npf, nfe, ne);
        }
      }
                  
      PGEMNMStridedBached(handle, m12, m, n, one, res.Ki, m12, DinvF, n, one, res.Hi, m12, ne, backend);

      schurVectorRh(tmp.tempn, res.Ri, npf, ncu12, ne); 
      ArrayCopy(res.Ri, tmp.tempn, m12*ne);
      PGEMNMStridedBached(handle, m12, 1, n, minusone, res.Ki, m12, Ru, n, one, res.Ri, m12, ne, backend);                 
    }
    
    if (common.debugMode==1) {    
      string filename;
      if (common.mpiProcs==1)
        filename = common.fileout;
      else
        filename = common.fileout + NumberToString(common.mpiRank);      
      writearray2file(filename + "uEquationElemSchur_DinvF.bin", DinvF, n*m*ne, backend);
      writearray2file(filename + "uEquationElemSchur_Ru.bin", Ru, n*ne, backend);
      writearray2file(filename + "uEquationElemSchur_DinvH.bin", DinvH, m*m*ne, backend);  
      writearray2file(filename + "uEquationElemSchur_Rh.bin", Rh, m*ne, backend);      
    } 
}

template <typename Model>
void uEquationHDG(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common,         
        cublasHandle_t handle, Int backend)
{    
    for (Int j=0; j<common.nbe; j++) {         
        uEquationElemBlock<Model>(sol, res, app, master, mesh, tmp, common, handle, j, backend);
        uEquationElemFaceBlock<Model>(sol, res, app, master, mesh, tmp, common, handle, j, backend);
        uEquationSchurBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
    }                     
}

template <typename Model>
void RuEquationElemBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int jth, Int backend)
{        
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    //Int ncq = common.ncq;// number of compoments of (q)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg) 
    Int ncs = common.ncs;// number of compoments of (sdg) 
    Int ncw = common.ncw;// number of compoments of (wdg) 
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe; // number of nodes on master element
    Int nge = common.nge; // number of gauss points on master element        

    Int e1 = common.eblks[3*jth]-1;
    Int e2 = common.eblks[3*jth+1];            
    Int ne = e2-e1;
    Int nn =  npe*ne; 
    Int nga = nge*ne;   
    //Int n0 = 0;                      // xg
    Int n1 = nga*ncx;                  // Xx
    Int n2 = nga*(ncx+nd*nd);          // jac        
    Int n3 = 0;                        // fg    
    Int n4 = nga*ncu*nd;               // sg  
    Int n5 = n4 + nga*ncu;             // ug
    Int n6 = n5 + nga*nc;              // wg
    Int nm = nge*e1*(ncx+nd*nd+1);

    dstype *xg = &sol.elemg[nm];
    dstype *Xx = &sol.elemg[nm+n1];
    dstype *jac = &sol.elemg[nm+n2];
    dstype *og = &sol.odgg[nge*nco*e1];

    dstype *wsrcg = &tmp.tempg[n3];        
    dstype *fg = &tmp.tempg[n3];        
    dstype *sg = &tmp.tempg[n4];
    dstype *uqg = &tmp.tempg[n5];
    dstype *wg = &tmp.tempg[n6];                

    GetArrayAtIndex(tmp.tempn, sol.udg, &mesh.eindudg1[npe*nc*e1], nn*nc);
    Node2Gauss(handle, uqg, tmp.tempn, master.shapegt, nge, npe, ne*nc, backend);
    
    if ((ncw>0) & (common.wave==0)) {
        GetElemNodes(tmp.tempn, sol.wdg, npe, ncw, 0, ncw, e1, e2);    
        Node2Gauss(handle, wg, tmp.tempn, master.shapegt, nge, npe, ne*ncw, backend);        
        
        GetElemNodes(tmp.tempn, sol.wsrc, npe, ncw, 0, ncw, e1, e2);    
        Node2Gauss(handle, wsrcg, tmp.tempn, master.shapegt, nge, npe, ne*ncw, backend);        

        // solve the w equation to get wg 
        wEquation(wg, xg, uqg, og, wsrcg, tmp.tempn, app, common, nga, backend); // fix bug here        
//         print2darray(wg, 1, 10);
//         print2darray(uqg, 1, 10);        
//         //exp(w) - sym(1.0) - u*u
//         error("here");
    }
    
    if (common.tdep) { // for time-dependent problem                
        // calculate sdg = sdg-udg*dtfactor
        ArrayAXPBY(sg, &sol.sdgg[nge*ncs*e1], uqg, one, -common.dtfactor, nga*ncu);            
        
        if (common.tdfunc==1) {
            // calculate the time derivative function Tdfunc(xdg, udg, odg)
            TdfuncDriver<Model>(fg, xg, uqg, og, wg, mesh, master, app, sol, tmp, common, nge, e1, e2, backend);

            // calculate (sdg-udg*dtfactor)*Tdfunc(xdg, udg, odg) 
            ArrayAXY(sg, sg, fg, one, nga*ncu);                
        }
        
        if (common.source==1) {            
            ArraySetValue(fg, 0.0, nga*ncu);
            // calculate the source term Source(xdg, udg, odg, wdg)
            SourceDriver<Model>(fg, xg, uqg, og, wg, mesh, master, app, sol, tmp, common, nge, e1, e2, backend);

            // calculate Source(xdg, udg, odg) + (sdg-udg*dtfactor)*Tdfunc(xdg, udg, odg) 
            ArrayAXPBY(sg, sg, fg, one, one, nga*ncu);            
        }        
    }
    else {  // steady state problem      
        ArraySetValue(sg, 0.0, nga*ncu);
        // calculate the source term Source(xdg, udg, odg, wdg)
        SourceDriver<Model>(sg, xg, uqg, og, wg, mesh, master, app, sol, tmp, common, nge, e1, e2, backend);                 
    }                
                
    ArraySetValue(fg, 0.0, nga*ncu*nd);    
    FluxDriver<Model>(fg, xg, uqg, og, wg, mesh, master, app, sol, tmp, common, nge, e1, e2, backend);    
    
    // Ru = npe * ne * ncu 
    ApplyXxJac(tmp.tempn, sg, fg, Xx, jac, nge, nd, ncu, ne);
    Gauss2Node(handle, &res.Ru[npe*ncu*e1], tmp.tempn, master.shapegw, nge*(nd+1), npe, ncu*ne, backend); // fixed bug here                   
}

template <typename Model>
void RuEquationElemFaceBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int jth, Int backend)
{            
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    //Int ncq = common.ncq;// number of compoments of (q)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int ncw = common.ncw;
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe; // number of nodes on master element
    Int npf = common.npf; // number of nodes on master face           
    Int ngf = common.ngf; // number of gauss poInts on master face              
    Int nfe = common.nfe; // number of faces in each element

    Int e1 = common.eblks[3*jth]-1;
    Int e2 = common.eblks[3*jth+1];            
    Int ne = e2-e1;
    Int nf = nfe*ne;
    Int nn =  npf*nf; 
    Int nga = ngf*nf;   
    Int n0 = 0;                                 // xg
    Int n1 = nga*ncx;                           // nlg
    Int n2 = nga*(ncx+nd);                      // jac
    //Int n3 = nga*(0);                           // uhg    
    Int n4 = nga*(ncu);                         // udg
    Int n5 = nga*(ncu+nc);                      // odg
    Int n6 = nga*(ncu+nc+nco);                  // wsrc
    Int n7 = nga*(ncu+nc+nco+ncw);              // wdg
    Int n8 = nga*(ncu+nc+nco+ncw+ncw);          // fhg
    Int nm = ngf*nfe*e1*(ncx+nd+1);
    
    dstype *xg = &sol.elemfaceg[nm+n0];    
    dstype *nlg = &sol.elemfaceg[nm+n1];
    dstype *jac = &sol.elemfaceg[nm+n2];

    dstype *uhg = &tmp.tempg[0];
    dstype *udg = &tmp.tempg[n4];
    dstype *odg = &tmp.tempg[n5];    
    dstype *wdg = &tmp.tempg[n6];
    dstype *wsrcg = &tmp.tempg[n7];
    dstype *fh     = &tmp.tempg[n8];
    
    // uhg = tmp.tempg[n3] at gauss points on face
    GetElementFaceNodes(tmp.tempn, sol.uh, mesh.elemcon, npf*nfe, ncu, e1, e2, 0); // fixed bug here

    // udg = tmp.tempg[n4] at gauss points on face
    GetElementFaceNodes(&tmp.tempn[nn*ncu], sol.udg, mesh.perm, npf*nfe, nc, npe, nc, e1, e2);

    if (nco>0) GetElementFaceNodes(&tmp.tempn[nn*(ncu+nc)], sol.odg, mesh.perm, npf*nfe, nco, npe, nco, e1, e2);      

    // fix bug here
    if ((ncw>0) & (common.wave==0)) {
      GetElementFaceNodes(&tmp.tempn[nn*(ncu+nc+nco)], sol.wdg, mesh.perm, npf*nfe, ncw, npe, ncw, e1, e2); 
      GetElementFaceNodes(&tmp.tempn[nn*(ncu+nc+nco+ncw)], sol.wsrc, mesh.perm, npf*nfe, ncw, npe, ncw, e1, e2); 
    }
    
    Node2Gauss(handle, tmp.tempg, tmp.tempn, master.shapfgt, ngf, npf, nfe*ne*(ncu+nc+nco+ncw+ncw), backend);
    
    if ((ncw>0) & (common.wave==0)) {
        // copy udg to tmp.tempn
        ArrayCopy(tmp.tempn, udg, nga*nc);
        
        // replace u with uhat 
        ArrayCopy(tmp.tempn, uhg, nga*ncu);
          
        // solve the w equation to get wg 
        wEquation(wdg, xg, tmp.tempn, odg, wsrcg, &tmp.tempn[nga*nc], app, common, nga, backend);                
        
        // solve the w equation to get wg 
        // wEquation(wdg, xg, udg, odg, wsrcg, tmp.tempn, app, common, nga, backend);                
    }
    
    FhatDriver<Model>(fh, tmp.tempn, xg, udg, odg, wdg, uhg, nlg, mesh, master, app, sol, tmp, common, nga, backend);      
    columnwiseMultiply(fh, fh, jac, nga, ncu);

    dstype *Rutmp = &tmp.tempn[0];
    Gauss2Node(handle, Rutmp, fh, master.shapfgw, ngf, npf, nf*ncu, backend);                
    assembleRu(&res.Ru[npe*ncu*e1], Rutmp, mesh.perm, npe, npf*nfe, ne*ncu);

    // impose boundary conditions
    for (int ibc=0; ibc<common.maxnbc; ibc++)
    {
      int n = ibc + common.maxnbc*jth;
      int start = common.nboufaces[n];
      int nfaces = common.nboufaces[n + 1] - start;
      if (nfaces>0) {
        int ngb = nfaces*ngf;
        dstype *xgb = &tmp.tempg[n8];
        dstype *ugb = &tmp.tempg[n8 + ngb*ncx];
        dstype *ogb = &tmp.tempg[n8 + ngb*ncx + ngb*nc];
        dstype *wgb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco]; 
        dstype *uhb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw];
        dstype *nlb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu];
        dstype *wsb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu + ngb*nd];
        dstype *fhb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu + ngb*nd + ngb*ncw];
        dstype *Rb =  &tmp.tempn[npf*nfe*ne*ncu];

        GetBoundaryNodes(xgb, xg, &mesh.boufaces[start], ngf, nfe, ne, ncx, nfaces);
        GetBoundaryNodes(ugb, udg, &mesh.boufaces[start], ngf, nfe, ne, nc, nfaces);
        GetBoundaryNodes(ogb, odg, &mesh.boufaces[start], ngf, nfe, ne, nco, nfaces);
        GetBoundaryNodes(wgb, wdg, &mesh.boufaces[start], ngf, nfe, ne, ncw, nfaces);
        GetBoundaryNodes(wsb, wsrcg, &mesh.boufaces[start], ngf, nfe, ne, ncw, nfaces);
        GetBoundaryNodes(uhb, uhg, &mesh.boufaces[start], ngf, nfe, ne, ncu, nfaces);
        GetBoundaryNodes(nlb, nlg, &mesh.boufaces[start], ngf, nfe, ne, nd, nfaces);

        if ((ncw>0) & (common.wave==0)) {
          // copy ugb to tmp.tempn
          ArrayCopy(Rb, ugb, ngb*nc);
        
          // replace u with uhat 
          ArrayCopy(Rb, uhb, ngb*ncu);
          
          wEquation(wgb, xgb, Rb, ogb, wsb, &Rb[ngb*nc], app, common, ngb, backend);          
          //wEquation(wgb, xgb, ugb, ogb, wsb, Rb, app, common, ngb, backend);
        }
        
        if (ibc+1 == 1000) StgInflowHDG(fhb, &tmp.tempg[n8], xgb, ogb, uhb, app.physicsparam, app.stgdata, 
                                 app.stgparam, common.time, ngb, common.stgNmode, nd);          
        else
            FbouDriver<Model>(fhb, xgb, ugb, ogb, wgb, uhb, nlb, 
             mesh, master, app, sol, tmp, common, ngb, ibc+1, backend);    

        dstype *jacb = &tmp.tempg[n8];
        GetBoundaryNodes(jacb, jac, &mesh.boufaces[start], ngf, nfe, ne, 1, nfaces);
        columnwiseMultiply(fhb, fhb, jacb, ngb, ncu);

        Gauss2Node(handle, Rb, fhb, master.shapfgw, ngf, npf, nfaces*ncu, backend);                
        PutBoundaryNodes(Rutmp, Rb, &mesh.boufaces[start], npf, nfe, ne, ncu, nfaces);        
      }
    }

    // npf*nfe*ne*ncu -> npf*nfe*ne*ncu
    ArrayAXPB(&res.Rh[npf*nfe*ncu*e1], Rutmp, minusone, zero, npf*nfe*ne*ncu);

    // permute Ru from npe*ne*ncu to npe*ncu*ne
    schurVectorRu(tmp.tempn, &res.Ru[npe*ncu*e1], npe, ncu, ne);
    ArrayCopy(&res.Ru[npe*ncu*e1], tmp.tempn, npe*ncu*ne);

    // permute Rh from npf*nfe*ne*ncu to ncu*npf*nfe*ne
    schurVectorRh(tmp.tempn, &res.Rh[npf*nfe*ncu*e1], npf*nfe, ncu, ne);
    ArrayCopy(&res.Rh[npf*nfe*ncu*e1], tmp.tempn, npf*nfe*ncu*ne);
    
    // impose interface conditions
    for (int ibc=0; ibc<common.maxnbc; ibc++)
    {
      int n = ibc + common.maxnbc*jth;
      int start = common.nboufaces[n];
      int nfaces = common.nboufaces[n + 1] - start;
      if ((nfaces>0) && (common.eblks[3*jth+2]==-1) && (common.coupledboundarycondition == ibc+1) && (common.coupledcondition>0))
      {
        int ncu12 = common.szinterfacefluxmap;
        int ngb = nfaces*ngf;
        dstype *xgb = &tmp.tempg[n8];
        dstype *ugb = &tmp.tempg[n8 + ngb*ncx];
        dstype *ogb = &tmp.tempg[n8 + ngb*ncx + ngb*nc];
        dstype *wgb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco]; 
        dstype *uhb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw];
        dstype *nlb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu];
        dstype *wsb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu + ngb*nd];
        dstype *fhb = &tmp.tempg[n8 + ngb*ncx + ngb*nc + ngb*nco + ngb*ncw + ngb*ncu + ngb*nd + ngb*ncw];
        dstype *Rb =  &tmp.tempn[npf*nfe*ne*ncu];

        GetBoundaryNodes(xgb, xg, &mesh.boufaces[start], ngf, nfe, ne, ncx, nfaces);
        GetBoundaryNodes(ugb, udg, &mesh.boufaces[start], ngf, nfe, ne, nc, nfaces);
        GetBoundaryNodes(ogb, odg, &mesh.boufaces[start], ngf, nfe, ne, nco, nfaces);
        GetBoundaryNodes(wgb, wdg, &mesh.boufaces[start], ngf, nfe, ne, ncw, nfaces);
        GetBoundaryNodes(wsb, wsrcg, &mesh.boufaces[start], ngf, nfe, ne, ncw, nfaces);
        GetBoundaryNodes(uhb, uhg, &mesh.boufaces[start], ngf, nfe, ne, ncu, nfaces);
        GetBoundaryNodes(nlb, nlg, &mesh.boufaces[start], ngf, nfe, ne, nd, nfaces);

        if ((ncw>0) & (common.wave==0)) {
          // copy ugb to tmp.tempn
          ArrayCopy(Rb, ugb, ngb*nc);
        
          // replace u with uhat 
          ArrayCopy(Rb, uhb, ngb*ncu);
          
          wEquation(wgb, xgb, &tmp.tempn[npf*nfe*ne*ncu], ogb, wsb, &Rb[ngb*nc], app, common, ngb, backend);          
          //wEquation(wgb, xgb, ugb, ogb, wsb, Rb, app, common, ngb, backend);
        }
        
        FintDriver<Model>(fhb, xgb, ugb, ogb, wgb, uhb, nlb, 
             mesh, master, app, sol, tmp, common, ngb, common.coupledcondition, backend);    

        dstype *jacb = &tmp.tempg[n8];
        GetBoundaryNodes(jacb, jac, &mesh.boufaces[start], ngf, nfe, ne, 1, nfaces);
        columnwiseMultiply(fhb, fhb, jacb, ngb, ncu12);

        Gauss2Node(handle, Rb, fhb, master.shapfgw, ngf, npf, nfaces*ncu12, backend);                
        ArrayAXPB(res.Ri, Rb, minusone, zero, npf*nfaces*ncu12);       
                
      }
    }            
}

template <typename Model>
void ResidualHDG(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common,         
        cublasHandle_t handle, Int backend)
{    
    for (Int j=0; j<common.nbe; j++) {        
        RuEquationElemBlock<Model>(sol, res, app, master, mesh, tmp, common, handle, j, backend);
        RuEquationElemFaceBlock<Model>(sol, res, app, master, mesh, tmp, common, handle, j, backend);        
    }                     
}


#endif

