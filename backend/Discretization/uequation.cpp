#ifndef __UEQUATION
#define __UEQUATION

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
    Int nn =  npe*ne; 
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
        wEquation(wg, wg_uq, xg, uqg, og, wsrc, app, common, tmp, nga, backend);
    }
    
    if (common.tdep) { // for time-dependent problem                
        // calculate sdg = sdg-udg*dtfactor
        ArrayAXPBY(sg, &sol.sdgg[nge*ncs*e1], uqg, one, -common.dtfactor, nga*ncu);            
        
        if (common.tdfunc==1) {
            // calculate the time derivative function Tdfunc(xdg, udg, odg)
            TdfuncDriver(fg, xg, uqg, og, wg, mesh, master, app, sol, tmp, common, nge, e1, e2, backend);

            // calculate (sdg-udg*dtfactor)*Tdfunc(xdg, udg, odg) 
            ArrayAXY(sg, sg, fg, one, nga*ncu);                
        }
        
        if (common.source==1) {            
            // calculate the source term Source(xdg, udg, odg, wdg)
            ArraySetValue(fg, 0.0, nga*ncu);
            ArraySetValue(sg_uq, 0.0, nga*ncu*nc);
            if ((ncw>0) & (common.wave==0)) ArraySetValue(sg_w, 0.0, nga*ncu*ncw);
            SourceDriver(fg, sg_uq, sg_w, xg, uqg, og, wg, mesh, master, app, sol, tmp, common, nge, e1, e2, backend);

            // calculate Source(xdg, udg, odg) + (sdg-udg*dtfactor)*Tdfunc(xdg, udg, odg) 
            ArrayAXPBY(sg, sg, fg, one, one, nga*ncu);            
        }        
    }
    else {  // steady state problem      
        // calculate the source term Source(xdg, udg, odg, wdg)
        ArraySetValue(sg, 0.0, nga*ncu);
        ArraySetValue(sg_uq, 0.0, nga*ncu*nc);
        if ((ncw>0) & (common.wave==0)) ArraySetValue(sg_w, 0.0, nga*ncu*ncw);
        SourceDriver(sg, sg_uq, sg_w, xg, uqg, og, wg, mesh, master, app, sol, tmp, common, nge, e1, e2, backend);                 
    }                
        
    // fg = nga*ncu*nd, sg = nga*ncu, fg_uq = nga*ncu*nd*nc, sg_uq = nga*ncu*nc

    ArraySetValue(fg, 0.0, nga*ncu*nd);
    ArraySetValue(fg_uq, 0.0, nga*ncu*nd*nc);
    if ((ncw>0) & (common.wave==0)) ArraySetValue(fg_w, 0.0, nga*ncu*nd*ncw); 
    FluxDriver(fg, fg_uq, fg_w, xg, uqg, og, wg, mesh, master, app, sol, tmp, common, nge, e1, e2, backend);    
    
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

    // dstype nrm[10];    
    // nrm[0] = PNORM(common.cublasHandle, npe*ncu*ne, &res.Ru[npe*ncu*e1], backend);       
    // nrm[1] = PNORM(common.cublasHandle, npe*npe*ncu*ncq*ne, res.B, backend);            
    // nrm[2] = PNORM(common.cublasHandle, npe*npe*ncu*ncu*ne, res.D, backend);            
    // cout<<"npe = "<<npe<<", ncu = "<<ncu<<", ncq = "<<ncq<<", ne = "<<ne<<endl;
    // cout<<"||Ru|| = "<<nrm[0]<<", ||B|| = "<<nrm[1]<<", ||D|| = "<<nrm[2]<<endl; 

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
    Int n3 = nga*(0);                           // uhg    
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
    }
    
    Node2Gauss(handle, tmp.tempg, tmp.tempn, master.shapfgt, ngf, npf, nfe*ne*(ncu+nc+nco+ncw), backend);

    if ((ncw>0) & (common.wave==0)) {
        // solve the w equation to get wg and wg_uq
        wEquation(wdg, wdg_uq, xg, udg, odg, wsrc, app, common, tmp, nga, backend);
    }

    ArraySetValue(fh, 0.0, nga*ncu);
    ArraySetValue(fh_uq, 0.0, nga*ncu*nc);
    ArraySetValue(fh_uh, 0.0, nga*ncu*ncu);
    if (ncw > 0) ArraySetValue(fh_w, 0.0, nga*ncu*ncw);       
    FhatDriver(fh, fh_uq, fh_w, fh_uh, xg, udg, odg, wdg, uhg, nlg, 
        mesh, master, app, sol, tmp, common, nga, backend);      

    if ((ncw>0) & (common.wave==0)) {
      ArrayGemmBatch2(fh_uq, fh_w, wdg_uq, one, ncu*nd, nc, ncw, nga);        
    }

    // if (common.debugMode==1) {    
    //   writearray2file(common.fileout + "uEquationElemFace_fh.bin", fh, ngf*nfe*ne*ncu, backend);
    // }

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
      //if (common.mpiRank==0) cout<<ibc<<" "<<n<<" "<<start<<" "<<nfaces<<endl;
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
        dstype *Rb =  &tmp.tempn[npf*nfe*ne*ncu + npf*npf*nfe*ne*ncu*ncu + npf*npf*nfe*ne*ncu*ncq + npf*npf*nfe*ne*ncu*ncu];

        GetBoundaryNodes(xgb, xg, &mesh.boufaces[start], ngf, nfe, ne, ncx, nfaces);
        GetBoundaryNodes(ugb, udg, &mesh.boufaces[start], ngf, nfe, ne, nc, nfaces);
        GetBoundaryNodes(ogb, odg, &mesh.boufaces[start], ngf, nfe, ne, nco, nfaces);
        GetBoundaryNodes(wsb, wsrc, &mesh.boufaces[start], ngf, nfe, ne, ncw, nfaces);
        GetBoundaryNodes(uhb, uhg, &mesh.boufaces[start], ngf, nfe, ne, ncu, nfaces);
        GetBoundaryNodes(nlb, nlg, &mesh.boufaces[start], ngf, nfe, ne, nd, nfaces);

        if ((ncw>0) & (common.wave==0)) {
          wEquation(wgb, wgb_uq, xgb, ugb, ogb, wsb, app, common, tmp, ngb, backend);
        }
        
        // intialize fhb, fhb_uq, fhb_w, fhb_uh to zero 
        ArraySetValue(fhb, 0.0, ngb*ncu);
        ArraySetValue(fhb_uq, 0.0, ngb*ncu*nc);
        ArraySetValue(fhb_uh, 0.0, ngb*ncu*ncu);
        if (ncw > 0) ArraySetValue(fhb_w, 0.0, ngb*ncu*ncw);        
        FbouDriver(fhb, fhb_uq, fhb_w, fhb_uh, xgb, ugb, ogb, wgb, uhb, nlb, 
             mesh, master, app, sol, tmp, common, ngb, ibc+1, backend);    

        //print2darray(xgb, ngb, ncx);     
        // print2darray(nlb, ngb, nd);     
        // print2darray(ugb, ngb, nc);     
        // print2darray(uhb, ngb, ncu);     
        // print2darray(fhb_uq, ngb*ncu, nc);     

        if ((ncw>0) & (common.wave==0)) {
          ArrayGemmBatch2(fhb_uq, fhb_w, wgb_uq, one, ncu*nd, nc, ncw, ngb);        
        }

        dstype *jacb = &tmp.tempg[n8];
        GetBoundaryNodes(jacb, jac, &mesh.boufaces[start], ngf, nfe, ne, 1, nfaces);
        columnwiseMultiply(fhb, fhb, jacb, ngb, ncu);
        columnwiseMultiply(fhb_uq, fhb_uq, jacb, ngb, ncu*nc);
        columnwiseMultiply(fhb_uh, fhb_uh, jacb, ngb, ncu*ncu);

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

    // dstype nrm[10];    
    // nrm[0] = PNORM(common.cublasHandle, npe*ncu*ne, &res.Ru[npe*ncu*e1], backend);       
    // nrm[1] = PNORM(common.cublasHandle, npe*npe*ncu*ncq*ne, res.B, backend);            
    // nrm[2] = PNORM(common.cublasHandle, npe*npe*ncu*ncu*ne, res.D, backend);            
    // nrm[3] = PNORM(common.cublasHandle, npf*nfe*ne*ncu, &res.Rh[npf*nfe*ncu*e1], backend);       
    // nrm[4] = PNORM(common.cublasHandle,  npe*npf*nfe*ncu*ncu*ne, &res.F[npe*npf*nfe*ncu*ncu*e1], backend);            
    // nrm[5] = PNORM(common.cublasHandle,  npf*nfe*npe*ne*ncu*ncu, res.K, backend);            
    // nrm[6] = PNORM(common.cublasHandle,  npf*nfe*npe*ne*ncu*ncq, res.G, backend);            
    // nrm[7] = PNORM(common.cublasHandle,  npf*nfe*npf*nfe*ne*ncu*ncu, &res.H[npf*nfe*npf*nfe*ncu*ncu*e1], backend);            

    // cout<<"e1 = "<<e1<<", e2 = "<<e2<<endl;
    // cout<<"npe = "<<npe<<", npf = "<<npf<<", nfe = "<<nfe<<", ncu = "<<ncu<<", ncq = "<<ncq<<", ne = "<<ne<<endl;
    // cout<<"||Ru|| = "<<nrm[0]<<", ||B|| = "<<nrm[1]<<", ||D|| = "<<nrm[2]<<endl; 
    // cout<<"||Rh|| = "<<nrm[3]<<", ||F|| = "<<nrm[4]<<", ||K|| = "<<nrm[5]<<endl; 
    // cout<<"||G|| = "<<nrm[6]<<", ||H|| = "<<nrm[7]<<endl; 

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
    dstype *DinvH =&res.H[m*m*e1];
    dstype *Rh = &res.Rh[m*e1];

    dstype *D = tmp.tempn;
    dstype *F = &tmp.tempn[npe*ncu*npe*ncu*ne];
    dstype *K = &tmp.tempn[npe*ncu*npe*ncu*ne + npe*npf*nfe*ncu*ncu*ne];
    dstype *H = &tmp.tempn[npe*ncu*npe*ncu*ne + npe*npf*nfe*ncu*ncu*ne + npe*npf*nfe*ncu*ncu*ne];

    // npe * npe * ne * ncu * ncu -> npe * ncu * npe * ncu * ne
    schurMatrixD(D, res.D, npe, ncu, ne);
    // npe * npf * nfe * ne * ncu * ncu -> npe * ncu * ncu * npf * nfe * ne
    schurMatrixF(F, &res.F[n*m*e1], npe, ncu, npf, nfe, ne);
    // npf * nfe * npe * ne * ncu * ncu -> ncu * npf * nfe * npe * ncu * ne
    schurMatrixK(K, res.K, npe, ncu, npf, nfe, ne);    
    // npf * nfe * npf * nfe * ne * ncu * ncu -> ncu * npf * nfe * ncu * npf * nfe * ne
    schurMatrixH(H, &res.H[m*m*e1], ncu, npf, nfe, ne);

    // dstype nrm[10];    
    // nrm[0] = PNORM(common.cublasHandle, npe*npe*ncu*ncu*ne, D, backend);       
    // nrm[1] = PNORM(common.cublasHandle, npe*npf*nfe*ncu*ncu*ne, F, backend);       
    // nrm[2] = PNORM(common.cublasHandle,  npe*npf*nfe*ne*ncu*ncu, K, backend);            
    // nrm[3] = PNORM(common.cublasHandle,  npf*nfe*npf*nfe*ne*ncu*ncu, H, backend);            

    // cout<<"||D|| = "<<nrm[0]<<", ||F|| = "<<nrm[1]<<", ||K|| = "<<nrm[2]<<", ||H|| = "<<nrm[3]<<endl;     

    // nrm[0] = PNORM(common.cublasHandle, npe*npe*common.ne, res.C, backend);       
    // nrm[1] = PNORM(common.cublasHandle, npe*npe*common.ne, &res.C[npe*npe*common.ne], backend);       
    // nrm[2] = PNORM(common.cublasHandle,  npe*npf*nfe*common.ne, res.E, backend);            
    // nrm[3] = PNORM(common.cublasHandle,  npe*npf*nfe*common.ne, &res.E[npe*npf*nfe*common.ne], backend);            

    // cout<<"||Cx|| = "<<nrm[0]<<", ||Cy|| = "<<nrm[1]<<", ||Ex|| = "<<nrm[2]<<", ||Ey|| = "<<nrm[3]<<endl;     

    if (common.ncq > 0) {      
      if (nd == 1) {
        // D = D + B * Minv * C
        schurMatrixBMinvC(D, res.B, &res.C[npe*npe*e1], npe, ncu, ne);
        // F = F - B * Minv * E
        schurMatrixBMinvE(F, res.B, &res.E[npe*npf*nfe*e1], npe, ncu, npf, nfe, ne);
        // K = K + G * Minv * C
        schurMatrixGMinvC(K, res.G, &res.C[npe*npe*e1], npe, ncu, npf, nfe, ne);
        // H = H - G * Minv * E
        schurMatrixGMinvE(H, res.G, &res.E[npe*npf*nfe*e1], npe, ncu, npf, nfe, ne);
      } 
      else if (nd == 2) {
        dstype *Cx = &res.C[npe*npe*e1]; // fixed bug here
        dstype *Cy = &res.C[npe*npe*common.ne + npe*npe*e1]; // fixed bug here
        dstype *Ex = &res.E[npe*npf*nfe*e1]; // fixed bug here
        dstype *Ey = &res.E[npe*npf*nfe*common.ne + npe*npf*nfe*e1]; // fixed bug here
        dstype *Bx = res.B;
        dstype *By = &res.B[npe*npe*ncu*ncu*ne];
        dstype *Gx = res.G;
        dstype *Gy = &res.G[npf*nfe*npe*ncu*ncu*ne];

        // nrm[0] = PNORM(common.cublasHandle, npe*npe*ne, Cx, backend);       
        // nrm[1] = PNORM(common.cublasHandle, npe*npe*ne, Cy, backend);       
        // nrm[2] = PNORM(common.cublasHandle,  npe*npf*nfe*ne, Ex, backend);            
        // nrm[3] = PNORM(common.cublasHandle,  npe*npf*nfe*ne, Ey, backend);            

        // cout<<"||Cx|| = "<<nrm[0]<<", ||Cy|| = "<<nrm[1]<<", ||Ex|| = "<<nrm[2]<<", ||Ey|| = "<<nrm[3]<<endl;     

        schurMatrixBMinvC(D, Bx, Cx, npe, ncu, ne);
        schurMatrixBMinvC(D, By, Cy, npe, ncu, ne);

        schurMatrixBMinvE(F, Bx, Ex, npe, ncu, npf, nfe, ne);
        schurMatrixBMinvE(F, By, Ey, npe, ncu, npf, nfe, ne);

        schurMatrixGMinvC(K, Gx, Cx, npe, ncu, npf, nfe, ne);
        schurMatrixGMinvC(K, Gy, Cy, npe, ncu, npf, nfe, ne);

        schurMatrixGMinvE(H, Gx, Ex, npe, ncu, npf, nfe, ne); 
        schurMatrixGMinvE(H, Gy, Ey, npe, ncu, npf, nfe, ne);
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

        schurMatrixBMinvC(D, Bx, Cx, npe, ncu, ne);
        schurMatrixBMinvC(D, By, Cy, npe, ncu, ne);
        schurMatrixBMinvC(D, Bz, Cz, npe, ncu, ne);

        schurMatrixBMinvE(F, Bx, Ex, npe, ncu, npf, nfe, ne);
        schurMatrixBMinvE(F, By, Ey, npe, ncu, npf, nfe, ne);
        schurMatrixBMinvE(F, Bz, Ez, npe, ncu, npf, nfe, ne);

        schurMatrixGMinvC(K, Gx, Cx, npe, ncu, npf, nfe, ne);
        schurMatrixGMinvC(K, Gy, Cy, npe, ncu, npf, nfe, ne);
        schurMatrixGMinvC(K, Gz, Cz, npe, ncu, npf, nfe, ne);

        schurMatrixGMinvE(H, Gx, Ex, npe, ncu, npf, nfe, ne);
        schurMatrixGMinvE(H, Gy, Ey, npe, ncu, npf, nfe, ne);
        schurMatrixGMinvE(H, Gz, Ez, npe, ncu, npf, nfe, ne);
      }
    }

    // nrm[0] = PNORM(common.cublasHandle, npe*npe*ncu*ncu*ne, D, backend);       
    // nrm[1] = PNORM(common.cublasHandle, npe*npf*nfe*ncu*ncu*ne, F, backend);       
    // nrm[2] = PNORM(common.cublasHandle,  npe*npf*nfe*ne*ncu*ncu, K, backend);            
    // nrm[3] = PNORM(common.cublasHandle,  npf*nfe*npf*nfe*ne*ncu*ncu, H, backend);            

    // cout<<"||D|| = "<<nrm[0]<<", ||F|| = "<<nrm[1]<<", ||K|| = "<<nrm[2]<<", ||H|| = "<<nrm[3]<<endl;     

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
    Inverse(handle, D, res.D, res.ipiv, n, ne, backend);    

    // DinvF = -Dinv * F
    PGEMNMStridedBached(handle, n, m, n, minusone, D, n, F, n, zero, DinvF, n, ne, backend); // fixed bug here

    // Ru = Dinv * Ru    
    schurVectorRu(F, Ru, npe, ncu, ne); // permute Ru from npe*ne*ncu to npe*ncu*ne
    PGEMNMStridedBached(handle, n, 1, n, one, D, n, F, n, zero, Ru, n, ne, backend);    

    // H = H + K * DinvF
    ArrayCopy(DinvH, H, m*m*ne);
    PGEMNMStridedBached(handle, m, m, n, one, K, m, DinvF, n, one, DinvH, m, ne, backend);

    // Rh = Rh - K * Dinv * Ru
    schurVectorRh(F, Rh, npf*nfe, ncu, ne); // permute Rh from npf*nfe*ne*ncu to ncu*npf*nfe*ne
    ArrayCopy(Rh, F, m*ne);
    PGEMNMStridedBached(handle, m, 1, n, minusone, K, m, Ru, n, one, Rh, m, ne, backend);     

    // nrm[0] = PNORM(common.cublasHandle, npe*ncu*ne, Ru, backend);       
    // nrm[1] = PNORM(common.cublasHandle, npf*nfe*ne*ncu, Rh, backend);       
    // nrm[2] = PNORM(common.cublasHandle,  npe*npf*nfe*ncu*ncu*ne, DinvF, backend);            
    // nrm[3] = PNORM(common.cublasHandle,  npf*nfe*npf*nfe*ne*ncu*ncu, DinvH, backend);            
    // cout<<"||Ru|| = "<<nrm[0]<<", ||Rh|| = "<<nrm[1]<<", ||F|| = "<<nrm[2]<<", ||H|| = "<<nrm[3]<<endl;     

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

void uEquationHDG(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common,         
        cublasHandle_t handle, Int backend)
{    
    for (Int j=0; j<common.nbe; j++) {                
        uEquationElemBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
        uEquationElemFaceBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
        uEquationSchurBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
    }                     
}

void RuEquationElemBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
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

    dstype *wsrc = &tmp.tempg[n3];        
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
        Node2Gauss(handle, wsrc, tmp.tempn, master.shapegt, nge, npe, ne*ncw, backend);        

        // solve the w equation to get wg 
        wEquation(wg, xg, uqg, og, wsrc, app, common, tmp, nga, backend);
    }
    
    if (common.tdep) { // for time-dependent problem                
        // calculate sdg = sdg-udg*dtfactor
        ArrayAXPBY(sg, &sol.sdgg[nge*ncs*e1], uqg, one, -common.dtfactor, nga*ncu);            
        
        if (common.tdfunc==1) {
            // calculate the time derivative function Tdfunc(xdg, udg, odg)
            TdfuncDriver(fg, xg, uqg, og, wg, mesh, master, app, sol, tmp, common, nge, e1, e2, backend);

            // calculate (sdg-udg*dtfactor)*Tdfunc(xdg, udg, odg) 
            ArrayAXY(sg, sg, fg, one, nga*ncu);                
        }
        
        if (common.source==1) {            
            ArraySetValue(fg, 0.0, nga*ncu);
            // calculate the source term Source(xdg, udg, odg, wdg)
            SourceDriver(fg, xg, uqg, og, wg, mesh, master, app, sol, tmp, common, nge, e1, e2, backend);

            // calculate Source(xdg, udg, odg) + (sdg-udg*dtfactor)*Tdfunc(xdg, udg, odg) 
            ArrayAXPBY(sg, sg, fg, one, one, nga*ncu);            
        }        
    }
    else {  // steady state problem      
        ArraySetValue(sg, 0.0, nga*ncu);
        // calculate the source term Source(xdg, udg, odg, wdg)
        SourceDriver(sg, xg, uqg, og, wg, mesh, master, app, sol, tmp, common, nge, e1, e2, backend);                 
    }                
                
    ArraySetValue(fg, 0.0, nga*ncu*nd);    
    FluxDriver(fg, xg, uqg, og, wg, mesh, master, app, sol, tmp, common, nge, e1, e2, backend);    
    
    // Ru = npe * ne * ncu 
    ApplyXxJac(tmp.tempn, sg, fg, Xx, jac, nge, nd, ncu, ne);
    Gauss2Node(handle, &res.Ru[npe*ncu*e1], tmp.tempn, master.shapegw, nge*(nd+1), npe, ncu*ne, backend); // fixed bug here                   
}

void RuEquationElemFaceBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
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
    Int n3 = nga*(0);                           // uhg    
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
    
    // uhg = tmp.tempg[n3] at gauss points on face
    GetElementFaceNodes(tmp.tempn, sol.uh, mesh.elemcon, npf*nfe, ncu, e1, e2, 0); // fixed bug here

    // udg = tmp.tempg[n4] at gauss points on face
    GetElementFaceNodes(&tmp.tempn[nn*ncu], sol.udg, mesh.perm, npf*nfe, nc, npe, nc, e1, e2);

    if (nco>0) GetElementFaceNodes(&tmp.tempn[nn*(ncu+nc)], sol.odg, mesh.perm, npf*nfe, nco, npe, nco, e1, e2);      

    if ((ncw>0) & (common.wave==0)) {
      GetElementFaceNodes(&tmp.tempn[nn*(ncu+nc+nco)], sol.wsrc, mesh.perm, npf*nfe, ncw, npe, ncw, e1, e2); 
    }
    
    Node2Gauss(handle, tmp.tempg, tmp.tempn, master.shapfgt, ngf, npf, nfe*ne*(ncu+nc+nco+ncw), backend);

    if ((ncw>0) & (common.wave==0)) {
        // solve the w equation to get wg 
        wEquation(wdg, xg, udg, odg, wsrc, app, common, tmp, nga, backend);
    }
    
    FhatDriver(fh, tmp.tempn, xg, udg, odg, wdg, uhg, nlg, mesh, master, app, sol, tmp, common, nga, backend);      
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
        dstype *Rb =  &tmp.tempn[npf*nfe*ne*ncu + npf*npf*nfe*ne*ncu*ncu + npf*npf*nfe*ne*ncu*ncq + npf*npf*nfe*ne*ncu*ncu];

        GetBoundaryNodes(xgb, xg, &mesh.boufaces[start], ngf, nfe, ne, ncx, nfaces);
        GetBoundaryNodes(ugb, udg, &mesh.boufaces[start], ngf, nfe, ne, nc, nfaces);
        GetBoundaryNodes(ogb, odg, &mesh.boufaces[start], ngf, nfe, ne, nco, nfaces);
        GetBoundaryNodes(wsb, wsrc, &mesh.boufaces[start], ngf, nfe, ne, ncw, nfaces);
        GetBoundaryNodes(uhb, uhg, &mesh.boufaces[start], ngf, nfe, ne, ncu, nfaces);
        GetBoundaryNodes(nlb, nlg, &mesh.boufaces[start], ngf, nfe, ne, nd, nfaces);

        if ((ncw>0) & (common.wave==0)) {
          wEquation(wgb, xgb, ugb, ogb, wsb, app, common, tmp, ngb, backend);
        }
        
        FbouDriver(fhb, xgb, ugb, ogb, wgb, uhb, nlb, 
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
    schurVectorRu(tmp.tempg, &res.Ru[npe*ncu*e1], npe, ncu, ne);
    ArrayCopy(&res.Ru[npe*ncu*e1], tmp.tempg, npe*ncu*ne);

    // permute Rh from npf*nfe*ne*ncu to ncu*npf*nfe*ne
    schurVectorRh(tmp.tempg, &res.Rh[npf*nfe*ncu*e1], npf*nfe, ncu, ne);
    ArrayCopy(&res.Rh[npf*nfe*ncu*e1], tmp.tempg, npf*nfe*ncu*ne);
}

void ResidualHDG(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common,         
        cublasHandle_t handle, Int backend)
{    
    for (Int j=0; j<common.nbe; j++) {        
        RuEquationElemBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
        RuEquationElemFaceBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);        
    }                     
}

// void uEquationFaceBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
//         meshstruct &mesh, tempstruct &tmp, commonstruct &common, 
//         cublasHandle_t handle, Int f1, Int f2, Int ib, Int backend)
// {            
//     Int nc = common.nc; // number of compoments of (u, q, p)
//     Int ncu = common.ncu;// number of compoments of (u)
//     Int nco = common.nco;// number of compoments of (o)
//     Int ncx = common.ncx;// number of compoments of (xdg)        
//     Int ncw = common.ncw;
//     Int nd = common.nd;     // spatial dimension    
//     Int npe = common.npe; // number of nodes on master element
//     Int npf = common.npf; // number of nodes on master face           
//     Int ngf = common.ngf; // number of gauss poInts on master face              

//     Int nf = f2-f1;
//     Int nn =  npf*nf; 
//     Int nga = ngf*nf;   
//     Int n0 = 0;                                 // xg
//     Int n1 = nga*ncx;                           // nlg
//     Int n2 = nga*(ncx+nd);                      // jac
//     Int n3 = nga*(0);                           // uhg    
//     Int n4 = nga*(ncu);                         // ug1
//     Int n5 = nga*(ncu+nc);                      // wg1
//     Int n6 = nga*(ncu+nc+ncw);                  // ug2
//     Int n7 = nga*(ncu+nc+ncw+nc);               // wg2
//     Int n8 = nga*(ncu+2*nc+2*ncw);              // fhg
//     //Int n7 = nga*(ncu+2*nc+ncw);                // wdg
//     Int nm = ngf*f1*(ncx+nd+1);
    
//     dstype *uhg = &tmp.tempg[0];
//     dstype *ug1 = &tmp.tempg[n4];
//     dstype *wg1 = &tmp.tempg[n5];
//     dstype *ug2 = &tmp.tempg[n6];
//     dstype *wg2 = &tmp.tempg[n7];
//     dstype *og1 = &sol.og1[ngf*nco*f1];
//     dstype *og2 = &sol.og2[ngf*nco*f1];
//     dstype *xg = &sol.faceg[nm+n0];    
//     dstype *nlg = &sol.faceg[nm+n1];

//     dstype *fh1 = &tmp.tempn[0];
//     dstype *fh2 = &tmp.tempn[nga*ncu];
//     dstype *fh1_uq = &tmp.tempn[2*nga*ncu];
//     dstype *fh2_uq = &tmp.tempn[2*nga*ncu + nga*ncu*nc];
//     dstype *fh1_uh = &tmp.tempn[2*nga*ncu + 2*nga*ncu*nc];
//     dstype *fh2_uh = &tmp.tempn[2*nga*ncu + 2*nga*ncu*nc + nga*ncu*ncu];
//     dstype *fh1_w =  &tmp.tempn[2*nga*ncu + 2*nga*ncu*nc + 2*nga*ncu*ncu];
//     dstype *fh2_w =  &tmp.tempn[2*nga*ncu + 2*nga*ncu*nc + 2*nga*ncu*ncu + nga*ncu*ncw];
    
//     // uhg = tmp.tempg[n3] at gauss points on face
//     GetElemNodes(tmp.tempn, sol.uh, npf, ncu, 0, ncu, f1, f2);
    
//     // ug1 = tmp.tempg[n4] at gauss points on face
//     //GetFaceNodes(tmp.tempn, sol.udg, mesh.facecon, npf, nc, npe, nc, f1, f2, 1, backend);      
//     GetArrayAtIndex(&tmp.tempn[nn*ncu], sol.udg, &mesh.findudg1[npf*nc*f1], nn*nc);    
//     if (ncw>0)
//         GetFaceNodes(&tmp.tempn[nn*(ncu+nc)], sol.wdg, mesh.facecon, npf, ncw, npe, ncw, f1, f2, 1);      
    
//     if (ib==0) {
//         // ug2 = tmp.tempg[n6] at gauss points on face
//         //GetFaceNodes(tmp.tempn, sol.udg, mesh.facecon, npf, nc, npe, nc, f1, f2, 2, backend);     
//         GetArrayAtIndex(&tmp.tempn[nn*(ncu+nc+ncw)], sol.udg, &mesh.findudg2[npf*nc*f1], nn*nc);
//         if (ncw>0)
//             GetFaceNodes(&tmp.tempn[nn*(ncu+nc+ncw+nc)], sol.wdg, mesh.facecon, npf, ncw, npe, ncw, f1, f2, 2);        
//         Node2Gauss(handle, tmp.tempg, tmp.tempn, master.shapfgt, ngf, npf, nf*(ncu+2*nc+2*ncw), backend);
//     }
//     else
//         Node2Gauss(handle, tmp.tempg, tmp.tempn, master.shapfgt, ngf, npf, nf*(ncu+nc+ncw), backend);
        
//     // calculate fhat
//     if (ib==0) { // interior faces                
//         FhatDriver(fh1, fh1_uq, fh1_w, fh1_uh, xg, ug1, og1, wg1, uhg, nlg, 
//             mesh, master, app, sol, tmp, common, nga, backend);      

//         // form BD, F, Ru, GK, H, and Rh           

//         FhatDriver(fh2, fh2_uq, fh2_w, fh2_uh, xg, ug2, og2, wg2, uhg, nlg, 
//             mesh, master, app, sol, tmp, common, nga, backend);        

//         // form BD, F, Ru, GK, H, and Rh       
//     }
//     else { // boundary faces      
//         FhatDriver(fh1, fh1_uq, fh1_w, fh1_uh, xg, ug1, og1, wg1, uhg, nlg, 
//             mesh, master, app, sol, tmp, common, nga, backend);      

//         // form BD, F, and Ru           

//         FbouDriver(fh2, fh2_uq, fh2_w, fh2_uh, xg, ug1, og1, wg1, uhg, nlg, 
//             mesh, master, app, sol, tmp, common, nga, ib, backend);    

//         // form GK, H, and Rh    
//     }        
            
//     // // evaluate fhg * jac at gauss points on face
//     // ApplyJacFhat(&tmp.tempg[n3], &tmp.tempg[n8], &sol.faceg[nm+n2], nga, ncu, ngf);    
    
//     // // <fhat, w>_F = <jac fhat, w>_T = w * (fhg * jac): npf*ncu*nf        
//     // Gauss2Node(handle, &res.Rh[npf*ncu*f1], &tmp.tempg[n3], master.shapfgw, ngf, npf, nf*ncu, backend);            
    
// #ifdef EXADEBUG                           
//     writearray2file(common.fileout + NumberToString(ib) + "RuFace_uhgf.bin", &tmp.tempg[n3], ngf*ncu*nf, backend);  
//     writearray2file(common.fileout + NumberToString(ib) + "RuFace_fgf.bin", &tmp.tempg[n4], ngf*ncu*nf, backend);  
//     writearray2file(common.fileout + NumberToString(ib) + "RuFace_rnf.bin", tmp.tempn, npf*ncu*nf, backend);
//     writearray2file(common.fileout + NumberToString(ib) + "RuFace_ruf.bin", res.Ruf, npe*ncu*common.ne1, backend);
// #endif              
// }

// void uEquationFace(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
//         meshstruct &mesh, tempstruct &tmp, commonstruct &common,
//         cublasHandle_t handle, Int nbf1, Int nbf2, Int backend)
// {    
//     for (Int j=nbf1; j<nbf2; j++) {
//         Int f1 = common.fblks[3*j]-1;
//         Int f2 = common.fblks[3*j+1];    
//         Int ib = common.fblks[3*j+2];    
//         uEquationFaceBlock(sol, res, app, master, mesh, tmp, common, handle, f1, f2, ib, backend);
//     }                          
// }

#endif

