#ifndef __QOICALCULATION
#define __QOICALCULATION

void qoiElemBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
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
    Int nga = nge*ne;   
    Int n1 = nga*ncx;                  // Xx
    Int n2 = nga*(ncx+nd*nd);          // jac        
    Int nm = nge*e1*(ncx+nd*nd+1);

    dstype *xg = &sol.elemg[nm];
    dstype *Xx = &sol.elemg[nm+n1];
    dstype *jac = &sol.elemg[nm+n2];

    dstype *og = &sol.odgg[nge*nco*e1];
    dstype *uqg = &tmp.tempg[0];
    dstype *wg = &tmp.tempg[nga*nc];    
    dstype *sg = &tmp.tempg[nga*(nc+ncw)];    
            
    GetElemNodes(tmp.tempn, sol.udg, npe, nc, 0, nc, e1, e2);   
    Node2Gauss(handle, uqg, tmp.tempn, master.shapegt, nge, npe, ne*nc, backend);        
    if ((ncw>0) & (common.wave==0)) {
        GetElemNodes(tmp.tempn, sol.wdg, npe, ncw, 0, ncw, e1, e2);    
        Node2Gauss(handle, wg, tmp.tempn, master.shapegt, nge, npe, ne*ncw, backend);        
    }
    
    int nvqoi = common.nvqoi;     
    ArraySetValue(sg, 0.0, nga*nvqoi);     
    QoIvolumeDriver(sg, xg, uqg, og, wg, mesh, master, app, sol, tmp, common, nge, e1, e2, backend);    
    
    ApplyJac(sg, jac, nge*ne, nge*ne*nvqoi);     
    Gauss2Node(handle, tmp.tempn, sg, master.gwe, nge, 1, nvqoi*ne, backend); 
    
    ArraySetValue(tmp.tempg, 1.0, ne);
    for (int i = 0; i<nvqoi; i++) {
        dstype dotprod = 0;
        PDOT(handle, ne, tmp.tempg, inc1, &tmp.tempn[i*ne], inc1, &dotprod, backend);
        common.qoivolume[i] += dotprod;
    }    
}

void qoiElement(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common)
{    
    for (int i = 0; i<common.nvqoi; i++) common.qoivolume[i] = 0.0;
    for (Int j=0; j<common.nbe; j++) {              
        Int e2 = common.eblks[3*j+1];            
        if (e2 <= common.ne1) qoiElemBlock(sol, res, app, master, mesh, tmp, common, common.cublasHandle, j, common.backend);        
    }                     
}

void qoiFaceBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
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
    Int nm = ngf*f1*(ncx+nd+1);
    Int n1 = nga*ncx;                           // nlg
    Int n2 = nga*(ncx+nd);                      // jac

    GetElemNodes(tmp.tempn, sol.uh, npf, ncu, 0, ncu, f1, f2);        
    //GetFaceNodes(tmp.tempn, sol.udg, mesh.facecon, npf, nc, npe, nc, f1, f2, 1, backend);      
    GetArrayAtIndex(&tmp.tempn[nn*ncu], sol.udg, &mesh.findudg1[npf*nc*f1], nn*nc);        
    if (ncw>0)
        GetFaceNodes(&tmp.tempn[nn*(ncu+nc)], sol.wdg, mesh.facecon, npf, ncw, npe, ncw, f1, f2, 1);             
    Node2Gauss(handle, tmp.tempg, tmp.tempn, master.shapfgt, ngf, npf, nf*(ncu+nc+ncw), backend);
        
    int nsurf = common.nsurf;     
    ArraySetValue(tmp.tempn, 0.0, nga*nsurf);     
    QoIboundaryDriver(tmp.tempn, &sol.faceg[nm], &tmp.tempg[nga*ncu], &sol.og1[ngf*nco*f1], 
            &tmp.tempg[nga*(ncu+nc)], &tmp.tempg[0], &sol.faceg[nm+n1], mesh, master, app, 
            sol, tmp, common, ngf, f1, f2, ib, backend);        

    ApplyJac(tmp.tempn, &sol.faceg[nm+n2], nga, nga*nsurf);
    Gauss2Node(handle, tmp.tempg, tmp.tempn, master.gwf, ngf, 1, nsurf*nf, backend); 

    ArraySetValue(tmp.tempn, 1.0, nf);
    for (int i = 0; i<nsurf; i++) {
        dstype dotprod = 0;
        PDOT(handle, nf, tmp.tempn, inc1, &tmp.tempg[i*nf], inc1, &dotprod, backend);
        common.qoisurface[i] += dotprod;
    }        
}

void qoiFace(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common)
{    
    for (int i = 0; i<common.nsurf; i++) common.qoisurface[i] = 0.0;
    for (Int j=0; j<common.nbf; j++) {
        Int f1 = common.fblks[3*j]-1;
        Int f2 = common.fblks[3*j+1];    
        Int ib = common.fblks[3*j+2];    
        if ((common.ibs > 0) && (ib == common.ibs))
            qoiFaceBlock(sol, res, app, master, mesh, tmp, common, common.cublasHandle, f1, f2, 1, common.backend);
    }                          
}

#endif

