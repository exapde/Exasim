#ifndef __GETUHAT
#define __GETUHAT

void UhatBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
        Int nd, Int npe, Int npf, Int nc, Int ncu, Int ncx, Int nco, Int f1, Int f2, 
        Int ib, Int jb, Int backend)
{        
    Int ncq = ncu*nd;
    Int nf = f2-f1;
    Int nn = npf*nf; 
    Int nga = npf*nf;   
    Int n0 = 0;                                 // xg
    Int n1 = nga*ncx;                           // nlg
    Int n2 = nga*(ncx+nd);                      // jac
    Int n3 = nga*(ncx+nd+1);                    // Jg, uhg
    Int n4 = nga*(ncx+nd+1+ncu);                // ug
    Int n5 = nga*(ncx+nd+1+ncu+nc);             // og
    
    if (ib==0) {        
        GetFaceNodes(tmp.tempn, sol.udg, mesh.facecon, npf, ncu, npe, nc, f1, f2, 0, backend);
        PutElemNodes(sol.uh, tmp.tempn, npf, ncu, 0, ncu, f1, f2, backend);
    }
    else {
        //GetFaceNodes(tmp.tempn, sol.xdg, mesh.facecon, npf, ncx, npe, ncx, f1, f2, 1, backend);            
        //Node2Gauss(handle, &tmp.tempg[n0], tmp.tempn, master.shapfnt, npf, npf, nf*ncx, backend);    
        GetArrayAtIndex(tmp.tempn, sol.xdg, &mesh.findxdg1[npf*ncx*f1], nn*ncx, backend);
        Node2Gauss(handle, &tmp.tempg[n0], tmp.tempn, master.shapfnt, npf, npf, nf*ncx, backend);    
    
        if (nd==1) {
            FaceGeom1D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga, backend);
        }
        else if (nd==2){
            Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapfnt[npf*npf], npf, npf, nf*nd, backend);                
            FaceGeom2D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga, backend);
        }
        else if (nd==3) {
            Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapfnt[npf*npf], npf, npf, nf*nd, backend);                     
            Node2Gauss(handle, &tmp.tempg[n3+nga*nd], tmp.tempn, &master.shapfnt[2*npf*npf], npf, npf, nf*nd, backend);                
            FaceGeom3D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga, backend);
        }        
        
        GetElemNodes(&tmp.tempg[n3], sol.uh, npf, ncu, 0, ncu, f1, f2, backend);
        //GetFaceNodes(&tmp.tempg[n4], sol.udg, mesh.facecon, npf, nc, npe, nc, f1, f2, 1, backend);
        GetArrayAtIndex(&tmp.tempg[n4], sol.udg, &mesh.findudg1[npf*nc*f1], nn*nc, backend);        
        //GetArrayAtIndex(tmp.tempn, sol.udg, &mesh.findudg1[npf*nc*f1], nn*nc, backend);
        //Node2Gauss(handle, &tmp.tempg[n4], tmp.tempn, master.shapfnt, npf, npf, nf*nc, backend);    
        
        //if (nco>0) {
        //    GetFaceNodes(&tmp.tempg[n5], sol.odg, mesh.facecon, npf, nco, npe, nco, f1, f2, 1,backend);       
        //}
        
        //cout<<common.tdepbcuh[jb]<<endl;
        if (common.tdep==0) {
            UbouDriver(tmp.tempn, &tmp.tempg[n0], &tmp.tempg[n3], &tmp.tempg[n4], &tmp.tempg[n5], 
                &tmp.tempg[n1], mesh, master, app, sol, tmp, common, npf, f1, f2, ib, backend);                    
        }
        else {
            if (common.tdepbcuh[jb]<0) 
                UbouDriver(tmp.tempn, &tmp.tempg[n0], &tmp.tempg[n3], &tmp.tempg[n4], &tmp.tempg[n5], 
                    &tmp.tempg[n1], mesh, master, app, sol, tmp, common, npf, f1, f2, ib, backend);                    
            else            
                UbouDriver2(tmp.tempn, &tmp.tempg[n0], &tmp.tempg[n3], &tmp.tempg[n4], &tmp.tempg[n5], 
                    &sol.udgbou0[common.tdepbcudg[jb]], &sol.udgbou1[common.tdepbcudg[jb]], &sol.udgbou2[common.tdepbcudg[jb]], 
                    &sol.uhbou0[common.tdepbcuh[jb]], &sol.uhbou1[common.tdepbcuh[jb]], &sol.uhbou2[common.tdepbcuh[jb]],
                    &tmp.tempg[n1], mesh, master, app, sol, tmp, common, npf, f1, f2, ib, jb, backend);            
        }
                                
        PutElemNodes(sol.uh, tmp.tempn, npf, ncu, 0, ncu, f1, f2, backend);                        
    }
// void UbouDriver2(dstype *fb, dstype *xg, dstype *uhg, dstype *udg, dstype * odg, dstype *nl, 
//         meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
//         commonstruct &common, Int ngf, Int f1, Int f2, Int ib, Int jb, Int backend)
    
}

// void UhatBlock2(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
//         meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
//         Int nd, Int npe, Int npf, Int nc, Int ncu, Int ncx, Int nco, Int f1, Int f2, Int ib, Int backend)
// {        
//     Int ncq = ncu*nd;
//     Int nf = f2-f1;
//     Int nn = npf*nf; 
//     Int nga = npf*nf;   
//     Int n0 = 0;                                 // xg
//     Int n1 = nga*ncx;                           // nlg
//     Int n2 = nga*(ncx+nd);                      // jac
//     Int n3 = nga*(ncx+nd+1);                    // Jg, uhg
//     Int n4 = nga*(ncx+nd+1+ncu);                // ug
//     Int n5 = nga*(ncx+nd+1+ncu+nc);             // og
//     Int nm = ngf*f1*(ncx+nd+1);
//     
//     if (ib==0) {        
//         GetFaceNodes(tmp.tempn, sol.udg, mesh.facecon, npf, ncu, npe, nc, f1, f2, 0, backend);
//         PutElemNodes(sol.uh, tmp.tempn, npf, ncu, 0, ncu, f1, f2, backend);
//     }
//     else {        
//         //GetElemNodes(&tmp.tempg[n3], sol.uh, npf, ncu, 0, ncu, f1, f2, backend);
//         
//         //GetFaceNodes(&tmp.tempg[n4], sol.udg, mesh.facecon, npf, nc, npe, nc, f1, f2, 1, backend);
//         GetArrayAtIndex(&tmp.tempg[n4], sol.udg, &mesh.findudg1[npf*nc*f1], nn*nc, backend);
//         
//         //if (nco>0) {
//         //    GetFaceNodes(&tmp.tempg[n5], sol.odg, mesh.facecon, npf, nco, npe, nco, f1, f2, 1,backend);       
//         //}
//         
//         UbouDriver(tmp.tempn, &sol.faceg[nm+n0], &tmp.tempg[n3], &tmp.tempg[n4], &tmp.tempg[n5], 
//                 &sol.faceg[nm+n1], mesh, master, app, sol, tmp, common, npf, f1, f2, ib, backend);
//                         
//         PutElemNodes(sol.uh, tmp.tempn, npf, ncu, 0, ncu, f1, f2, backend);                        
//     }
// }

void GetUhat(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
        Int nbf1, Int nbf2, Int backend)
{        
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)   
    Int nd = common.nd;     // spatial dimension        
    Int npe = common.npe; // number of nodes on master element
    Int npf = common.npf; // number of nodes on master face      
    
    for (Int j=nbf1; j<nbf2; j++) {
        Int f1 = common.fblks[3*j]-1;
        Int f2 = common.fblks[3*j+1];    
        Int ib = common.fblks[3*j+2];    
        //UhatBlock(sol, res, app, master, mesh, tmp, common, handle, f1, f2, ib, backend);
        UhatBlock(sol, res, app, master, mesh, tmp, common,
                handle, nd, npe, npf, nc, ncu, ncx, nco, f1, f2, ib, j, backend);
    }                           
}

void InitUhat(solstruct &sol, meshstruct &mesh, tempstruct &tmp, commonstruct &common, Int backend)
{        
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)   
    Int nd = common.nd;     // spatial dimension        
    Int npe = common.npe; // number of nodes on master element
    Int npf = common.npf; // number of nodes on master face      
    
    for (Int j=0; j<common.nbf; j++) {
        Int f1 = common.fblks[3*j]-1;
        Int f2 = common.fblks[3*j+1];    
        Int ib = common.fblks[3*j+2];    
        
        if (ib==0) {        
            GetFaceNodes(tmp.tempn, sol.udg, mesh.facecon, npf, ncu, npe, nc, f1, f2, 0, backend);
            PutElemNodes(sol.uh, tmp.tempn, npf, ncu, 0, ncu, f1, f2, backend);
        }        
        else {
            GetArrayAtIndex(tmp.tempn, sol.udg, &mesh.findudg1[npf*nc*f1], npf*(f2-f1)*nc, backend);
            PutElemNodes(sol.uh, tmp.tempn, npf, ncu, 0, ncu, f1, f2, backend);
        }        
    }                           
}

void SetUdgUhatBou(solstruct &sol, commonstruct &common, Int backend)
{        
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)   
    Int nd = common.nd;     // spatial dimension        
    Int npe = common.npe; // number of nodes on master element
    Int npf = common.npf; // number of nodes on master face      
    Int n1 = 0;
    Int n2 = 0;
        
    for (Int j=0; j<common.nbf; j++) {
        Int f1 = common.fblks[3*j]-1;
        Int f2 = common.fblks[3*j+1];    
        Int ib = common.fblks[3*j+2];    
        Int nf = f2-f1;
        Int nn = npf*nf;   
        common.tdepbcuh[j] = -1;
        common.tdepbcudg[j] = -1;
        for (Int k=0; k<common.ntdepbc; k++) {            
            if (ib == common.tdepbc[k]) {
                common.tdepbcuh[j] = n1;
                common.tdepbcudg[j] = n2;
                n1 = n1 + nn*ncu;
                n2 = n2 + nn*nc;
            }            
        }
    }   
    
    // allocate memory
    if (n1>0) {
        TemplateMalloc(&sol.uhbou0, n1, backend);       
        TemplateMalloc(&sol.udgbou0, n2, backend);       
        if (common.tstages>1) {
            TemplateMalloc(&sol.uhbou1, n1, backend);       
            TemplateMalloc(&sol.udgbou1, n2, backend);       
        }
        if (common.tstages>2) {
            TemplateMalloc(&sol.uhbou2, n1, backend);       
            TemplateMalloc(&sol.udgbou2, n2, backend);       
        }
    }    
}

void GetUdgUhatBou(dstype *uhbou, dstype *udgbou, solstruct &sol, meshstruct &mesh, commonstruct &common, Int backend)
{        
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)   
    Int nd = common.nd;     // spatial dimension        
    Int npe = common.npe; // number of nodes on master element
    Int npf = common.npf; // number of nodes on master face      
    
    for (Int j=0; j<common.nbf; j++) {
        Int f1 = common.fblks[3*j]-1;
        Int f2 = common.fblks[3*j+1];    
        Int ib = common.fblks[3*j+2];    
        Int nf = f2-f1;
        Int nn = npf*nf;           
        if (common.tdepbcuh[j]>=0) {
            Int n1 = common.tdepbcuh[j];
            Int n2 = common.tdepbcudg[j];                                
            GetElemNodes(&uhbou[n1], sol.uh, npf, ncu, 0, ncu, f1, f2, backend);
            GetArrayAtIndex(&udgbou[n2], sol.udg, &mesh.findudg1[npf*nc*f1], nn*nc, backend);                                            
        }
    }                           
}

#endif

