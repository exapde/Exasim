/*
    getuhat.cpp

    This file contains functions for computing the numerical trace (uhat) and its derivatives on the faces of elements
    in a finite element mesh, as part of the Exasim backend discretization routines.

    Functions:

    - int isin(Int ib, Int *a, Int n)
        Utility function to check if integer ib is present in array a of length n.

    - void UhatBlock(...)
        Computes the numerical trace (uhat) on a block of faces, handling both interior and boundary faces.
        For boundary faces, it computes geometric quantities, gathers solution and auxiliary data, and calls
        the boundary condition driver (UbouDriver). For interior faces, it directly copies solution data.

    - void GetUhat(...)
        Loops over face blocks and calls UhatBlock for each block to compute uhat for all faces in the specified range.

    - void dUhatBlock(...) [ifdef HAVE_ENZYME]
        Computes the derivative of the numerical trace (duhat) with respect to the solution variables, for use in
        automatic differentiation (Enzyme). Handles both interior and boundary faces, gathering necessary data and
        calling the boundary condition driver for derivatives.

    - void GetdUhat(...) [ifdef HAVE_ENZYME]
        Loops over face blocks and calls dUhatBlock for each block to compute duhat for all faces in the specified range.

    Notes:
    - The functions rely on various mesh, solution, and temporary data structures, as well as BLAS handles for
      possible GPU acceleration.
    - The code supports synthetic turbulence generation (commented out) and is designed for extensibility.
    - The HAVE_ENZYME preprocessor directive enables derivative computations for use with automatic differentiation.

    Parameters (common across functions):
    - solstruct, resstruct, appstruct, masterstruct, meshstruct, tempstruct, commonstruct: Data structures holding
      solution, mesh, application, master element, temporary, and common parameters.
    - cublasHandle_t handle: BLAS handle for GPU computations.
    - Int nd, npe, npf, nc, ncu, ncx, nco: Mesh and solution dimensions.
    - Int f1, f2, ib: Face block indices and boundary type.
    - Int backend: Backend identifier for CPU/GPU execution.

*/
#ifndef __GETUHAT
#define __GETUHAT

bool isin(Int ib, Int *a, Int n)
{
    bool in = false;
    for (int i=0; i<n; i++) {        
        if (ib == a[i]) {
            in = true;
            break;
        }
    }

    return in;
}        

template <typename Model>
void UhatBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
        Int nd, Int npe, Int npf, Int nc, Int ncu, Int ncx, Int nco, Int f1, Int f2, Int ib, Int backend)
{        
    Int ncw = common.ncw;
    //Int ncq = ncu*nd;
    Int nf = f2-f1;
    Int nn = npf*nf; 
    Int nga = npf*nf;   
    Int n0 = 0;                                 // xg
    Int n1 = nga*ncx;                           // nlg
    Int n2 = nga*(ncx+nd);                      // jac
    Int n3 = nga*(ncx+nd+1);                    // Jg, uhg
    Int n4 = nga*(ncx+nd+1+ncu);                // ug
    Int n5 = nga*(ncx+nd+1+ncu+nc);             // og
    Int n6 = nga*(ncx+nd+1+ncu+nc+ncw);         // wg
    
    if (ib==0) {        
        GetFaceNodes(tmp.tempn, sol.udg, mesh.facecon, npf, ncu, npe, nc, f1, f2, 0);
        PutElemNodes(sol.uh, tmp.tempn, npf, ncu, 0, ncu, f1, f2);
    }
    else if (isin(ib, common.stgib, common.nstgib)) {        
        //if (common.mpiRank==0 && ib > 0) printf("ib = %d \n", ib);
        dstype *xgb = &tmp.tempg[0];
        dstype *ogb = &tmp.tempg[nga*ncx];

        //GetFaceNodes(xgb, sol.xdg, mesh.facecon, npf, ncx, npe, ncx, f1, f2, 1);  
        GetArrayAtIndex(xgb, sol.xdg, &mesh.findxdg1[npf*ncx*f1], nn*ncx);
        if (nco>0) {
            GetFaceNodes(ogb, sol.odg, mesh.facecon, npf, ncu, npe, nco, f1, f2, 1);      
        }
        
        StgInflowLDG(tmp.tempn, xgb, ogb, app.physicsparam, app.stgdata, 
                          app.stgparam, common.time, nga, common.stgNmode, common.nd);          

        PutElemNodes(sol.uh, tmp.tempn, npf, ncu, 0, ncu, f1, f2);
    }
    //else if (isin(ib, common.stgib, common.nstgib)) {        
//     else if (ib <= 2) {        
//         // synthetic turbulence generation    
// //         GetArrayAtIndex(tmp.tempg, sol.xdg, &mesh.findxdg1[npf*ncx*f1], nn*ncx, backend);
// //         StgHomoTurb2(tmp.tempn, tmp.tempg, app.stgdata, app.uinf, app.stgparam, &app.stgparam[nd], 
// //                 app.physicsparam, common.time, nn, common.stgNmode, nd, backend);
//         
//         GetArrayAtIndex(tmp.tempn, sol.xdg, &mesh.findxdg1[npf*ncx*f1], nn*ncx, backend);
//         Node2Gauss(handle, &tmp.tempg[n0], tmp.tempn, master.shapfnt, npf, npf, nf*ncx, backend);    
//     
//         if (nd==1) {
//             FaceGeom1D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga, backend);
//         }
//         else if (nd==2){
//             Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapfnt[npf*npf], npf, npf, nf*nd, backend);                
//             FaceGeom2D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga, backend);
//         }
//         else if (nd==3) {
//             Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapfnt[npf*npf], npf, npf, nf*nd, backend);                     
//             Node2Gauss(handle, &tmp.tempg[n3+nga*nd], tmp.tempn, &master.shapfnt[2*npf*npf], npf, npf, nf*nd, backend);                
//             FaceGeom3D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga, backend);
//         }        
//         
//         // uh
//         StgHomoTurb(&tmp.tempg[n3], &tmp.tempg[n0], app.stgdata,  app.stgparam,
//                 common.time, nn, common.stgNmode, nd, backend);
// //         StgHomoTurb2(&tmp.tempg[n3], &tmp.tempg[n0], app.stgdata, app.uinf, app.stgparam, &app.stgparam[nd], 
// //                 app.physicsparam, common.time, nn, common.stgNmode, nd, backend);
//         
//         // udg
//         GetArrayAtIndex(&tmp.tempg[n4], sol.udg, &mesh.findudg1[npf*nc*f1], nn*nc, backend);        
//         
//         // odg
//         if (nco>0) {
//             GetFaceNodes(&tmp.tempg[n5], sol.odg, mesh.facecon, npf, nco, npe, nco, f1, f2, 1,backend);       
//         }
//         if (ncw>0) {
//             GetFaceNodes(&tmp.tempg[n6], sol.wdg, mesh.facecon, npf, ncw, npe, ncw, f1, f2, 1,backend);       
//         }
//         
//         UbouDriver<Model>(tmp.tempn, &tmp.tempg[n0], &tmp.tempg[n3], &tmp.tempg[n4], &tmp.tempg[n5], 
//                 &tmp.tempg[n6], &tmp.tempg[n1], mesh, master, app, sol, tmp, common, 
//                 npf, f1, f2, ib, backend);
//         
//         PutElemNodes(sol.uh, tmp.tempn, npf, ncu, 0, ncu, f1, f2, backend);                        
//         
// //         printArray2D(tmp.tempn,nn,ncu,backend);
// //         printArray2D(app.stgdata,common.stgNmode,10,backend);
// //         printArray2D(app.stgparam,1,2*nd,backend);
// //         printArray2D(app.uinf,1,4,backend);
// //         cout<<common.stgNmode<<endl;
// //         error("here");        
//     }
    else {        
        //GetFaceNodes(tmp.tempn, sol.xdg, mesh.facecon, npf, ncx, npe, ncx, f1, f2, 1, backend);            
        //Node2Gauss(handle, &tmp.tempg[n0], tmp.tempn, master.shapfnt, npf, npf, nf*ncx, backend);    
        GetArrayAtIndex(tmp.tempn, sol.xdg, &mesh.findxdg1[npf*ncx*f1], nn*ncx);
        Node2Gauss(handle, &tmp.tempg[n0], tmp.tempn, master.shapfnt, npf, npf, nf*ncx, backend);
    
        if (nd==1) {
            FaceGeom1D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga);
        }
        else if (nd==2){
            Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapfnt[npf*npf], npf, npf, nf*nd, backend);                
            FaceGeom2D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga);
        }
        else if (nd==3) {
            Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapfnt[npf*npf], npf, npf, nf*nd, backend);                     
            Node2Gauss(handle, &tmp.tempg[n3+nga*nd], tmp.tempn, &master.shapfnt[2*npf*npf], npf, npf, nf*nd, backend);                
            FaceGeom3D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga);
        }        
        
        //GetElemNodes(&tmp.tempg[n3], sol.uh, npf, ncu, 0, ncu, f1, f2, backend);
        //GetFaceNodes(&tmp.tempg[n4], sol.udg, mesh.facecon, npf, nc, npe, nc, f1, f2, 1, backend);
        GetArrayAtIndex(&tmp.tempg[n4], sol.udg, &mesh.findudg1[npf*nc*f1], nn*nc);
        
        if (nco>0) {
            GetFaceNodes(&tmp.tempg[n5], sol.odg, mesh.facecon, npf, nco, npe, nco, f1, f2, 1);       
        }
                
        UbouDriver<Model>(tmp.tempn, &tmp.tempg[n0], &tmp.tempg[n4], &tmp.tempg[n5], &tmp.tempg[n6], &tmp.tempg[n3], 
                 &tmp.tempg[n1], mesh, master, app, sol, tmp, common, npf, f1, f2, ib, backend);
                               
        PutElemNodes(sol.uh, tmp.tempn, npf, ncu, 0, ncu, f1, f2);
    }
}

template <typename Model>
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
        UhatBlock<Model>(sol, res, app, master, mesh, tmp, common,
                handle, nd, npe, npf, nc, ncu, ncx, nco, f1, f2, ib, backend);
    }                           
}

#ifdef HAVE_ENZYME
  //// Method 2
template <typename Model>
void dUhatBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
        Int nd, Int npe, Int npf, Int nc, Int ncu, Int ncx, Int nco, Int f1, Int f2, Int ib, Int backend)
{        
    Int ncw = common.ncw;
    //Int ncq = ncu*nd;
    Int nf = f2-f1;
    Int nn = npf*nf; 
    Int nga = npf*nf;   
    Int n0 = 0;                                 // xg
    Int n1 = nga*ncx;                           // nlg
    Int n2 = nga*(ncx+nd);                      // jac
    Int n3 = nga*(ncx+nd+1);                    // Jg, uhg
    Int n4 = nga*(ncx+nd+1+ncu);                // ug
    Int n5 = nga*(ncx+nd+1+ncu+nc);             // og
    Int n6 = nga*(ncx+nd+1+ncu+nc+ncw);         // wg
    Int n7 = nga*(ncx+nd+1+ncu+nc+ncw+ncw);     // duhg
    Int n8 = nga*(ncx+nd+1+ncu+nc+ncw+ncw+ncu); // dug
    Int n9 = nga*(ncx+nd+1+ncu+nc+ncw+ncw+ncu+nc); // dwg

    if (ib==0) { 
        GetFaceNodes(tmp.tempn, sol.udg, mesh.facecon, npf, ncu, npe, nc, f1, f2, 0);
        PutElemNodes(sol.uh, tmp.tempn, npf, ncu, 0, ncu, f1, f2);
        
        // calculates duh (duhat/du v) on interior
        // reuses tmp.tmpn array; v = sol.dudg
        GetFaceNodes(tmp.tempn, sol.dudg, mesh.facecon, npf, ncu, npe, nc, f1, f2, 0);
        PutElemNodes(sol.duh, tmp.tempn, npf, ncu, 0, ncu, f1, f2);
    } 
    else {
        // Geometry information
        GetArrayAtIndex(tmp.tempn, sol.xdg, &mesh.findxdg1[npf*ncx*f1], nn*ncx);
        Node2Gauss(handle, &tmp.tempg[n0], tmp.tempn, master.shapfnt, npf, npf, nf*ncx, backend);    
    
        if (nd==1) {
            FaceGeom1D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga);
        }
        else if (nd==2){
            Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapfnt[npf*npf], npf, npf, nf*nd, backend);                
            FaceGeom2D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga);
        }
        else if (nd==3) {
            Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapfnt[npf*npf], npf, npf, nf*nd, backend);                     
            Node2Gauss(handle, &tmp.tempg[n3+nga*nd], tmp.tempn, &master.shapfnt[2*npf*npf], npf, npf, nf*nd, backend);                
            FaceGeom3D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga);
        }        
        
        // udg = tmp.tmpg[n4] 
        GetArrayAtIndex(&tmp.tempg[n4], sol.udg, &mesh.findudg1[npf*nc*f1], nn*nc);
        // v = dudg = tmp.tmpg[n8]
        GetArrayAtIndex(&tmp.tempg[n8], sol.dudg, &mesh.findudg1[npf*nc*f1], nn*nc);

        // TODO: need to decide whether to set to 0 out here or inside drivers
        ArraySetValue(tmp.tempn, zero, 2*nn*ncu);

        // uhat = tempn[0], duhat = tempn[nn*ncu]
        UbouDriver<Model>(&tmp.tempn[0], &tmp.tempn[nn*ncu], &tmp.tempg[n0], &tmp.tempg[n4], &tmp.tempg[n8],
                &tmp.tempg[n5], &tmp.tempg[n6], &tmp.tempg[n9], &tmp.tempg[n3], 
                &tmp.tempg[n1], mesh, master, app, sol, tmp, common, npf, f1, f2, ib, backend);

        //Update sol.duh = (duhat/du) v   
        PutElemNodes(sol.duh, &tmp.tempn[nn*ncu], npf, ncu, 0, ncu, f1, f2);
    }
}

template <typename Model>
void GetdUhat(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
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
        dUhatBlock<Model>(sol, res, app, master, mesh, tmp, common,
                handle, nd, npe, npf, nc, ncu, ncx, nco, f1, f2, ib, backend);
    }                           
}


#endif

#endif

