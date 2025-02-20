#ifndef __WRESIDUAL
#define __WRESIDUAL

#include "../AppDriver/sourcewDriver.cpp"
#include "../AppDriver/eosDriver.cpp"
#include "../AppDriver/eosduDriver.cpp"
#include "../AppDriver/eosdwDriver.cpp"

// // Calculate Rwe = (u, nabla r)_K for a given u
// void RwElemBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
//         meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int nd, 
//         Int npe, Int nge, Int nc, Int ncu, Int ncx, Int e1, Int e2, Int backend)
// {            
//     
//     Int ncq = ncu*nd;
//     Int ne = e2-e1;
//     Int nn =  npe*ne; 
//     Int nga = nge*ne;   
//     Int n0 = 0;                                 // ug
//     Int n1 = nga*ncu;                           // fg
//     Int n2 = nge*e1*(ncx+nd*nd+1) + nga*ncx;
//         
//     // udg = tmp.tempg[n0] at gauss points on elements
//     GetArrayAtIndex(tmp.tempn, sol.udg, &mesh.eindudg1[npe*nc*e1], nn*nc, backend);
//     Node2Gauss(handle, &tmp.tempg[n0], tmp.tempn, master.shapegt, nge, npe, ne*ncu, backend);    
//     
//     // udg * Xx  at gauss points on elements: nge*nd*ncu*nd*ne   
//     ApplyXx3(&tmp.tempg[n1], &tmp.tempg[n0], &sol.elemg[n2], nge, nd, ncu, ne, backend);   
//     
//     // (u, nabla_i dot v)_K = (u, dv/dx_i)_K = (jac u, dv/dxi_j dxi_j/dx_i)_T = sum_j dv/dxi_j * udg * Xx(:,i,j)
//     Gauss2Node(handle, &res.Rqe[npe*ncq*e1], &tmp.tempg[n1], &master.shapegw[npe*nge], 
//             nge*nd, npe, ncq*ne, backend);                
// }
// 
// void RwElem(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
//         meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
//         Int nbe1, Int nbe2, Int backend)
// {
//     Int nc = common.nc; // number of compoments of (u, q, p)
//     Int ncu = common.ncu;// number of compoments of (u)
//     Int ncw = common.ncw;// number of compoments of (w)
//     Int nco = common.nco;// number of compoments of (o)
//     Int ncx = common.ncx;// number of compoments of (xdg)        
//     Int nd = common.nd;     // spatial dimension    
//     Int npe = common.npe; // number of nodes on master element
//     Int nge = common.nge; // number of gauss points on master element    
//     //Int ne = common.ne; // number of elements in this subdomain 
//     
//         
//     for (Int j=nbe1; j<nbe2; j++) {
//         Int e1 = common.eblks[3*j]-1;
//         Int e2 = common.eblks[3*j+1];    
//         if (common.dae==1)
//             SourcewDriver(&sol.wdg[npe*ncw*e1], &sol.xdg[npe*ncx*e1], &sol.udg[npe*nc*e1], &sol.odg[npe*nco*e1], 
//                     mesh, master, app, sol, tmp, common, npe, e1, e2, backend);
//         else
//             RwElemBlock(sol, res, app, master, mesh, tmp, common, common.cublasHandle, nd, npe, nge, nc, 
//                 ncu, ncx, e1, e2, backend);                
//     }                     
// }
// 
// // Calculate Rwf = <uhat dot n, r>_F for a given uhat
// void RwFaceBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
//         meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
//         Int nd, Int npe, Int npf, Int ngf, Int nc, Int ncu, Int ncx, Int f1, Int f2, Int ib, Int backend)
// {        
//     Int ncq = ncu*nd;
//     Int nf = f2-f1;
//     //Int nn =  npf*nf; 
//     Int nga = ngf*nf;   
//     Int n1 = nga*ncx;                           // nlg
//     Int n2 = nga*(ncx+nd);                      // jac
//     Int n3 = 0;                                 // uhg
//     Int n4 = nga*ncu;                           // fhg
// 
//     // uhg = tmp.tempg[n3] at gauss points on face
//     GetElemNodes(tmp.tempn, sol.uh, npf, ncu, 0, ncu, f1, f2, backend);      
//     Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, master.shapfgt, ngf, npf, nf*ncu, backend);
//     
//     // evaluate uhg * jac * nlg at gauss points on face
//     ApplyJacNormal(&tmp.tempg[n4], &tmp.tempg[n3], &sol.faceg[ngf*f1*(ncx+nd+1)+n1], 
//                 &sol.faceg[ngf*f1*(ncx+nd+1)+n2], nga, ncu, nd, ngf, backend);    
//     
//     // <uhat, v dot n>_F = <jac uhat, v dot n>_T = v * (uhg * jac * nlg)            
//     Gauss2Node(handle, &res.Rh[npf*ncq*f1], &tmp.tempg[n4], master.shapfgw, ngf, npf, nf*ncq, backend);            
//     
// #ifdef DEBUG                           
//     writearray2file(common.fileout + NumberToString(ib) + "RqFace_uhgf.bin", &tmp.tempg[n3], ngf*ncu*nf, backend);  
//     writearray2file(common.fileout + NumberToString(ib) + "RqFace_fgf.bin", &tmp.tempg[n4], ngf*ncq*nf, backend);  
//     writearray2file(common.fileout + NumberToString(ib) + "RqFace_rnf.bin", tmp.tempn, npf*ncq*nf, backend);
//     writearray2file(common.fileout + NumberToString(ib) + "RqFace_rqf.bin", res.Rqf, npe*ncq*common.ne1, backend);
// #endif              
// }
// 
// 
// void RwFace(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
//         meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int nbf1, Int nbf2, Int backend)
// {    
//     Int nc = common.nc; // number of compoments of (u, q, p)
//     Int ncu = common.ncu;// number of compoments of (u)
//     Int ncx = common.ncx;// number of compoments of (xdg)        
//     Int nd = common.nd;     // spatial dimension    
//     Int npe = common.npe; // number of nodes on master element
//     Int npf = common.npf; // number of nodes on master face           
//     Int ngf = common.ngf; // number of gauss poInts on master face          
//     //Int ne = common.ne; // number of elements in this subdomain 
//     
//     for (Int j=nbf1; j<nbf2; j++) {
//         Int f1 = common.fblks[3*j]-1;
//         Int f2 = common.fblks[3*j+1];    
//         Int ib = common.fblks[3*j+2];    
//         if (common.dae!=1)
//             RwFaceBlock(sol, res, app, master, mesh, tmp, common, common.cublasHandle, 
//                 nd, npe, npf, ngf, nc, ncu, ncx, f1, f2, ib, backend);
//     }                       
// }
// 

#endif

