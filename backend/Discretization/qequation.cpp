#ifndef __QEQUATION
#define __QEQUATION

void qEquationElem(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int backend)
{        
    Int ncx = common.ncx;// number of compoments of (xdg) 
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe; // number of nodes on master element
    Int nge = common.nge; // number of gauss points on master element        
    Int ne = common.ne; // number of elements in this subdomain
    Int nbe = common.nbe; // number of blocks for elements   
    Int neb = common.neb; // maximum number of elements per block

    TemplateMalloc(&res.Mass2, npe*npe*ne, backend);
    TemplateMalloc(&res.Minv2, npe*npe*ne, backend);
    TemplateMalloc(&res.C, npe*npe*ne*nd, backend);
    res.szC = npe*npe*ne*nd;
    res.szMinv2 = npe*npe*ne; 
    res.szMass2 = npe*npe*ne; 

    dstype *work=nullptr;  
    Int *ipiv=nullptr;
    TemplateMalloc(&work, npe*npe*neb, backend);      
    TemplateMalloc(&ipiv, npe+1, backend);           

//     print2darray(master.shapegw, npe, nge);
//     print2darray(master.shapegt, nge, npe);
//     print2darray(master.shapegwdotshapeg, npe, npe);

    for (Int j=0; j<nbe; j++) // for each block of elements
    {
        Int e1 = common.eblks[3*j]-1;
        Int e2 = common.eblks[3*j+1];
        Int ns = e2-e1;        
        Int nga = nge*ns;
        Int n1 = nga*ncx;                  // Xx
        Int n2 = nga*(ncx+nd*nd);          // jac        
        Int nm = nge*e1*(ncx+nd*nd+1);

        //dstype *xg = &sol.elemg[nm];
        dstype *Xx = &sol.elemg[nm+n1];
        dstype *jac = &sol.elemg[nm+n2];

        Gauss2Node(handle, &res.Minv2[npe*npe*e1], jac, master.shapegwdotshapeg, nge, npe*npe, ns, backend);                

        ArrayCopy(&res.Mass2[npe*npe*e1], &res.Minv2[npe*npe*e1], npe*npe*ns);

//         print2darray(jac, nge, ns);
//         cout<<e1<<" "<<e2<<endl;        
//         cout<<"Mass2: "<<endl;
//         print2darray(res.Mass2, npe, npe);
//         error("here");

        Inverse(handle, &res.Minv2[npe*npe*e1], work, ipiv, npe, ns, backend);

        if (nd==1) {
          Gauss2Node(handle, work, Xx, &master.shapegwdotshapeg[npe*npe*nge], nge, npe*npe, ns, backend);   
          PGEMNMStridedBached(handle, npe, npe, npe, one, &res.Minv2[npe*npe*e1], npe, work, npe, zero, &res.C[npe*npe*e1], npe, ns, backend); // fixed bug here
        }
        else if (nd==2) {          
          dstype *Cx = res.C;
          dstype *Cy = &res.C[npe*npe*ne];
          Gauss2Node(handle,  work, Xx, &master.shapegwdotshapeg[npe*npe*nge], nge, npe*npe, ns, backend);     
          Gauss2Node1(handle, work, &Xx[nga*2], &master.shapegwdotshapeg[npe*npe*nge*2], nge, npe*npe, ns, backend);     
          PGEMNMStridedBached(handle, npe, npe, npe, one, &res.Minv2[npe*npe*e1], npe, work, npe, zero, &Cx[npe*npe*e1], npe, ns, backend); // fixed bug here

          Gauss2Node(handle,  work, &Xx[nga*1], &master.shapegwdotshapeg[npe*npe*nge], nge, npe*npe, ns, backend);     
          Gauss2Node1(handle, work, &Xx[nga*3], &master.shapegwdotshapeg[npe*npe*nge*2], nge, npe*npe, ns, backend);  
          PGEMNMStridedBached(handle, npe, npe, npe, one, &res.Minv2[npe*npe*e1], npe, work, npe, zero, &Cy[npe*npe*e1], npe, ns, backend); // fixed bug here  
        }
        else if (nd==3) { 
          dstype *Cx = res.C;
          dstype *Cy = &res.C[npe*npe*ne];
          dstype *Cz = &res.C[npe*npe*ne*2];
          Gauss2Node(handle,  work, Xx, &master.shapegwdotshapeg[npe*npe*nge], nge, npe*npe, ns, backend);
          Gauss2Node1(handle, work, &Xx[nga*3], &master.shapegwdotshapeg[npe*npe*nge*2], nge, npe*npe, ns, backend);
          Gauss2Node1(handle, work, &Xx[nga*6], &master.shapegwdotshapeg[npe*npe*nge*3], nge, npe*npe, ns, backend);
          PGEMNMStridedBached(handle, npe, npe, npe, one, &res.Minv2[npe*npe*e1], npe, work, npe, zero, &Cx[npe*npe*e1], npe, ns, backend); // fixed bug here  

          Gauss2Node(handle,  work, &Xx[nga*1], &master.shapegwdotshapeg[npe*npe*nge], nge, npe*npe, ns, backend);
          Gauss2Node1(handle, work, &Xx[nga*4], &master.shapegwdotshapeg[npe*npe*nge*2], nge, npe*npe, ns, backend);
          Gauss2Node1(handle, work, &Xx[nga*7], &master.shapegwdotshapeg[npe*npe*nge*3], nge, npe*npe, ns, backend);
          PGEMNMStridedBached(handle, npe, npe, npe, one, &res.Minv2[npe*npe*e1], npe, work, npe, zero, &Cy[npe*npe*e1], npe, ns, backend); // fixed bug here    

          Gauss2Node(handle,  work, &Xx[nga*2], &master.shapegwdotshapeg[npe*npe*nge], nge, npe*npe, ns, backend);
          Gauss2Node1(handle, work, &Xx[nga*5], &master.shapegwdotshapeg[npe*npe*nge*2], nge, npe*npe, ns, backend);
          Gauss2Node1(handle, work, &Xx[nga*8], &master.shapegwdotshapeg[npe*npe*nge*3], nge, npe*npe, ns, backend);
          PGEMNMStridedBached(handle, npe, npe, npe, one, &res.Minv2[npe*npe*e1], npe, work, npe, zero, &Cz[npe*npe*e1], npe, ns, backend); // fixed bug here  
        }
    }

    if (common.debugMode==1) {

      string filename;
      if (common.mpiProcs==1)
        filename = common.fileout;
      else
        filename = common.fileout + NumberToString(common.mpiRank);      

      writearray2file(filename + "qEquationElem_Mass.bin", res.Mass2, npe*npe*ne, backend);  
      writearray2file(filename + "qEquationElem_Minv.bin", res.Minv2, npe*npe*ne, backend);
      writearray2file(filename + "qEquationElem_C.bin", res.C, npe*npe*ne*nd, backend);

    }

    TemplateFree(work, backend);
    TemplateFree(ipiv, backend);
}

// Calculate Rqf = <uhat, v dot n>_F for a given uhat
void qEquationFaceBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
        Int nd, Int nfe, Int npe, Int npf, Int ngf, Int ncx, Int f1, Int f2, Int ib, Int backend)
{        
    Int nf = f2-f1;
    Int nga = ngf*nf;   
    Int n1 = nga*ncx;                           // nlg
    Int n2 = nga*(ncx+nd);                      // jac

    dstype *nlg = &sol.faceg[ngf*f1*(ncx+nd+1)+n1];
    dstype *jac = &sol.faceg[ngf*f1*(ncx+nd+1)+n2];
    dstype *Etmp = &tmp.tempg[nga];
    
    if (nd==1) {
      ArrayAXY(tmp.tempg, nlg, jac, one, nga);      
      Gauss2Node(handle, Etmp, tmp.tempg, master.shapfgwdotshapfg, ngf, npf*npf, nf, backend);
      assembleMatrixE(res.E, Etmp, &mesh.facecon[2*npf*f1], &mesh.f2e[4*f1], npf, npe, nfe, f1, f2);
    }
    else if (nd==2) {
      dstype *Ex = res.E;
      dstype *Ey = &res.E[npe*npf*nfe*common.ne];

      ArrayAXY(tmp.tempg, nlg, jac, one, nga);
      Gauss2Node(handle, Etmp, tmp.tempg, master.shapfgwdotshapfg, ngf, npf*npf, nf, backend);
      assembleMatrixE(Ex, Etmp, &mesh.facecon[2*npf*f1], &mesh.f2e[4*f1], npf, npe, nfe, f1, f2);

      // print3darray(Etmp, npf, npf, nf);
      // print3darray(Etmp, npf, npf, nf);
      // print2darray(Ex, npe, npf*nfe);

      ArrayAXY(tmp.tempg, &nlg[nga], jac, one, nga);
      Gauss2Node(handle, Etmp, tmp.tempg, master.shapfgwdotshapfg, ngf, npf*npf, nf, backend);
      assembleMatrixE(Ey, Etmp, &mesh.facecon[2*npf*f1], &mesh.f2e[4*f1], npf, npe, nfe, f1, f2);
    }
    else if (nd==3) {
      dstype *Ex = res.E;
      dstype *Ey = &res.E[npe*npf*nfe*common.ne];
      dstype *Ez = &res.E[npe*npf*nfe*common.ne*2];

      ArrayAXY(tmp.tempg, nlg, jac, one, nga);
      Gauss2Node(handle, Etmp, tmp.tempg, master.shapfgwdotshapfg, ngf, npf*npf, nf, backend);
      assembleMatrixE(Ex, Etmp, &mesh.facecon[2*npf*f1], &mesh.f2e[4*f1], npf, npe, nfe, f1, f2);

      ArrayAXY(tmp.tempg, &nlg[nga], jac, one, nga);
      Gauss2Node(handle, Etmp, tmp.tempg, master.shapfgwdotshapfg, ngf, npf*npf, nf, backend);
      assembleMatrixE(Ey, Etmp, &mesh.facecon[2*npf*f1], &mesh.f2e[4*f1], npf, npe, nfe, f1, f2);

      ArrayAXY(tmp.tempg, &nlg[nga*2], jac, one, nga);
      Gauss2Node(handle, Etmp, tmp.tempg, master.shapfgwdotshapfg, ngf, npf*npf, nf, backend);
      assembleMatrixE(Ez, Etmp, &mesh.facecon[2*npf*f1], &mesh.f2e[4*f1], npf, npe, nfe, f1, f2);
    }        
}

void qEquationFace(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int backend)
{    
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe; // number of nodes on master element
    Int npf = common.npf; // number of nodes on master face           
    Int ngf = common.ngf; // number of gauss poInts on master face
    Int nbf = common.nbf; // number of blocks for faces 
    Int nfe = common.nfe; // number of faces per element
    Int ne = common.ne; // number of elements in this subdomain

    TemplateMalloc(&res.E, npe*npf*nfe*ne*nd, backend);         
    ArraySetValue(res.E, zero, npe*npf*nfe*ne*nd);
    res.szE = npe*npf*nfe*ne*nd; 
    
    for (Int j=0; j<nbf; j++) {
        Int f1 = common.fblks[3*j]-1;
        Int f2 = common.fblks[3*j+1];    
        Int ib = common.fblks[3*j+2];    
        cout<<f1<<" "<<f2<<" "<<ib<<endl;
        qEquationFaceBlock(sol, res, app, master, mesh, tmp, common, handle, 
                nd, nfe, npe, npf, ngf, ncx, f1, f2, ib, backend);
    }                       

    if (common.debugMode==1) {
      writearray2file(common.fileout + "qEquationFace_E.bin", res.E, npe*npf*nfe*nd*common.ne, backend);  
    }    
}

// Calculate Rqf = <uhat, v dot n>_F for a given uhat
void qEquationElemFaceBlock(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
        Int nd, Int nfe, Int npe, Int npf, Int ngf, Int ncx, Int e1, Int e2, Int backend)
{          
    Int ns = e2-e1;
    Int nga = ngf*nfe*ns;   
    Int n1 = nga*ncx;                           // nlg
    Int n2 = nga*(ncx+nd);                      // jac

    dstype *nlg = &sol.elemfaceg[ngf*nfe*e1*(ncx+nd+1)+n1];
    dstype *jac = &sol.elemfaceg[ngf*nfe*e1*(ncx+nd+1)+n2];
    //dstype *Etmp = &tmp.tempg[nga];
    
    dstype *work=nullptr;  
    dstype *Etmp=nullptr;
    TemplateMalloc(&work, npe*npf*nfe*ns, backend);      
    TemplateMalloc(&Etmp, npf*npf*nfe*ns, backend);           

    if (nd==1) {
      ArrayAXY(tmp.tempg, nlg, jac, one, nga);      
      Gauss2Node(handle, Etmp, tmp.tempg, master.shapfgwdotshapfg, ngf, npf*npf, nfe*ns, backend);
      ArraySetValue(work, zero, npe*npf*nfe*ns);
      assembleMatrixE(work, Etmp, mesh.perm, npf, npe, nfe, ns);
      PGEMNMStridedBached(handle, npe, npf*nfe, npe, one, &res.Minv2[npe*npe*e1], npe, work, npe, zero, &res.E[npe*npe*e1], npe, ns, backend); // fixed bug here  
//       printArray2D(work, npe, npf*nfe, backend);     
//       printArray2D(nlg, nga, ncx, backend);     
//       error("here");
    }
    else if (nd==2) {
      dstype *Ex = res.E;
      dstype *Ey = &res.E[npe*npf*nfe*common.ne];

      ArrayAXY(tmp.tempg, nlg, jac, one, nga);
      Gauss2Node(handle, Etmp, tmp.tempg, master.shapfgwdotshapfg, ngf, npf*npf, nfe*ns, backend);
      ArraySetValue(work, zero, npe*npf*nfe*ns);
      assembleMatrixE(work, Etmp, mesh.perm, npf, npe, nfe, ns);
      PGEMNMStridedBached(handle, npe, npf*nfe, npe, one, &res.Minv2[npe*npe*e1], npe, work, npe, zero, &Ex[npe*npf*nfe*e1], npe, ns, backend); // fixed bug here  

      // print2darray(jac, ngf, nfe);
      // print2darray(nlg, ngf, nfe);
      // print2darray(Etmp, npf, npf*nfe);      
      // print2darray(Ex, npe, npf*nfe);
      // print2iarray(mesh.perm, npf, nfe);

      ArrayAXY(tmp.tempg, &nlg[nga], jac, one, nga);
      Gauss2Node(handle, Etmp, tmp.tempg, master.shapfgwdotshapfg, ngf, npf*npf, nfe*ns, backend);
      ArraySetValue(work, zero, npe*npf*nfe*ns);
      assembleMatrixE(work, Etmp, mesh.perm, npf, npe, nfe, ns);
      PGEMNMStridedBached(handle, npe, npf*nfe, npe, one, &res.Minv2[npe*npe*e1], npe, work, npe, zero, &Ey[npe*npf*nfe*e1], npe, ns, backend); // fixed bug here  

      // print2darray(&nlg[nga], ngf, nfe);
      // print2darray(Etmp, npf, npf*nfe);      
      // print2darray(Ey, npe, npf*nfe);
    }
    else if (nd==3) {
      dstype *Ex = res.E;
      dstype *Ey = &res.E[npe*npf*nfe*common.ne];
      dstype *Ez = &res.E[npe*npf*nfe*common.ne*2];

      ArrayAXY(tmp.tempg, nlg, jac, one, nga);
      Gauss2Node(handle, Etmp, tmp.tempg, master.shapfgwdotshapfg, ngf, npf*npf, nfe*ns, backend);
      ArraySetValue(work, zero, npe*npf*nfe*ns);
      assembleMatrixE(work, Etmp, mesh.perm, npf, npe, nfe, ns);
      PGEMNMStridedBached(handle, npe, npf*nfe, npe, one, &res.Minv2[npe*npe*e1], npe, work, npe, zero, &Ex[npe*npf*nfe*e1], npe, ns, backend); // fixed bug here  
      
//       print2darray(jac, ngf, nfe);
//       print2darray(nlg, ngf, nfe);
//       print2darray(Etmp, npf, npf*nfe);      
//       print2darray(Ex, npe, npf*nfe);
//       print2iarray(mesh.perm, npf, nfe);
//       print2darray(master.shapfgwdotshapfg, npf*npf, ngf);
//       error("here");
      
      ArrayAXY(tmp.tempg, &nlg[nga], jac, one, nga);
      Gauss2Node(handle, Etmp, tmp.tempg, master.shapfgwdotshapfg, ngf, npf*npf, nfe*ns, backend);
      ArraySetValue(work, zero, npe*npf*nfe*ns);
      assembleMatrixE(work, Etmp, mesh.perm, npf, npe, nfe, ns);
      PGEMNMStridedBached(handle, npe, npf*nfe, npe, one, &res.Minv2[npe*npe*e1], npe, work, npe, zero, &Ey[npe*npf*nfe*e1], npe, ns, backend); // fixed bug here  

      ArrayAXY(tmp.tempg, &nlg[nga*2], jac, one, nga);
      Gauss2Node(handle, Etmp, tmp.tempg, master.shapfgwdotshapfg, ngf, npf*npf, nfe*ns, backend);
      ArraySetValue(work, zero, npe*npf*nfe*ns);
      assembleMatrixE(work, Etmp, mesh.perm, npf, npe, nfe, ns);
      PGEMNMStridedBached(handle, npe, npf*nfe, npe, one, &res.Minv2[npe*npe*e1], npe, work, npe, zero, &Ez[npe*npf*nfe*e1], npe, ns, backend); // fixed bug here  
      
//       print2darray(&Ex[npe*npf*nfe*e1], npe, npf*nfe);
//       print2darray(&Ey[npe*npf*nfe*e1], npe, npf*nfe);
//       print2darray(&Ez[npe*npf*nfe*e1], npe, npf*nfe);
//       cout<<e1<<endl;
//       error("here");
    }        

    TemplateFree(work, backend);
    TemplateFree(Etmp, backend);
}

void qEquationElemFace(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int backend)
{    
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe; // number of nodes on master element
    Int npf = common.npf; // number of nodes on master face           
    Int ngf = common.ngf; // number of gauss poInts on master face
    Int nbe = common.nbe; // number of blocks for elements 
    Int nfe = common.nfe; // number of faces per element
    Int ne = common.ne; // number of elements in this subdomain

    TemplateMalloc(&res.E, npe*npf*nfe*ne*nd, backend);         
    ArraySetValue(res.E, zero, npe*npf*nfe*ne*nd);
    res.szE = npe*npf*nfe*ne*nd; 
    
    for (Int j=0; j<nbe; j++) // for each block of elements
    {
        Int e1 = common.eblks[3*j]-1;
        Int e2 = common.eblks[3*j+1];
        qEquationElemFaceBlock(sol, res, app, master, mesh, tmp, common, handle, 
                nd, nfe, npe, npf, ngf, ncx, e1, e2, backend);
    }    

    if (common.debugMode==1) {
      string filename;
      if (common.mpiProcs==1)
        filename = common.fileout;
      else
        filename = common.fileout + NumberToString(common.mpiRank);      

      writearray2file(filename + "qEquationFace_E.bin", res.E, npe*npf*nfe*nd*common.ne, backend);  
    }    
}

void qEquation(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, Int backend)
{
    qEquationElem(sol, res, app, master, mesh, tmp, common, common.cublasHandle, backend);
    qEquationElemFace(sol, res, app, master, mesh, tmp, common, common.cublasHandle, backend);    
}

void hdgGetQ(dstype *udg, dstype *uhat, solstruct &sol, resstruct &res, meshstruct &mesh, tempstruct &tmp, commonstruct &common, Int backend)
{
    Int nc = common.nc;// number of compoments of (udg)        
    Int ncu = common.ncu; // number of compoments of (u)
    Int ncq = common.ncq; // number of compoments of (q)
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe; // number of nodes on master element
    Int npf = common.npf; // number of nodes on master face           
    //Int ngf = common.ngf; // number of gauss poInts on master face
    Int nbe = common.nbe; // number of blocks for elements 
    Int nfe = common.nfe; // number of faces per element
    Int ne = common.ne; // number of elements in this subdomain
    Int ndf = npf*nfe; // number of dofs on a face

    dstype scalar = 1.0;
    if (common.wave==1)
        scalar = 1.0/common.dtfactor;    
    
    // Solve: dtfactor * M * q = C * u - E * uhat + M * s
    //        q = (1/dtfactor) * (Minv * C * u - Minv * E * uhat + s)
    
    for (Int j=0; j<nbe; j++) // for each block of elements
    {
      Int e1 = common.eblks[3*j]-1;
      Int e2 = common.eblks[3*j+1];
      Int ns = e2 - e1;

      dstype *u = tmp.tempg;
      dstype *q = &tmp.tempg[npe*ncu*ns];           
      dstype *uh = tmp.tempn;
      ArrayExtract(u, udg, npe, nc, ne, 0, npe, 0, ncu, e1, e2);        // npe x ncu x ns
      GetElementFaceNodes(uh, uhat, mesh.elemcon, npf*nfe, ncu, e1, e2, 1); // npf*nfe x ncu x ns

      dstype *s = res.Rq;
      ArraySetValue(s, zero, npe*ncq*ns);
      
      // print3darray(uhat, common.ncu, common.npf, common.nf);
      //print3darray(uh, npf*nfe, common.ncu, ns);            
      dstype *Cx = res.C;        
      dstype *Ex = res.E;              
      if (nd==1) {        
        PGEMNMStridedBached(common.cublasHandle, npe, ncu, npe, one, &Cx[npe*npe*e1], npe, u, npe, zero, q, npe, ns, backend);
        PGEMNMStridedBached(common.cublasHandle, npe, ncu, npf*nfe, minusone, &Ex[npe*npf*nfe*e1], npe, uh, ndf, one, q, npe, ns, backend);        
        if (common.wave==1) {// get the source term due to the time derivative 
          ArrayExtract(s, sol.sdg, npe, nc, ne, 0, npe, ncu, ncu+ncq, e1, e2);          
          ArrayAXPBY(q, q, s, scalar, scalar, npe*ncu*ns); // scalar*(q + s)                   
        }
        ArrayInsert(udg, q, npe, nc, ne, 0, npe, ncu, 2*ncu, e1, e2);
      }
      else if (nd==2) {        
        dstype *Cy = &res.C[npe*npe*common.ne];
        dstype *Ey = &res.E[npe*npf*nfe*common.ne];
        
        PGEMNMStridedBached(common.cublasHandle, npe, ncu, npe, one, &Cx[npe*npe*e1], npe, u, npe, zero, q, npe, ns, backend);
        PGEMNMStridedBached(common.cublasHandle, npe, ncu, npf*nfe, minusone, &Ex[npe*npf*nfe*e1], npe, uh, ndf, one, q, npe, ns, backend);
        if (common.wave==1) {
          ArrayExtract(s, sol.sdg, npe, nc, ne, 0, npe, ncu, 2*ncu, e1, e2);          
          ArrayAXPBY(q, q, s, scalar, scalar, npe*ncu*ns); // scalar*(q + s)
        }
        ArrayInsert(udg, q, npe, nc, ne, 0, npe, ncu, 2*ncu, e1, e2);

        // print2darray(Cx, npe, npe);            
        // print2darray(Ex, npe, npf*nfe);            
        // print2darray(u, npe, ncu);            
        // print2darray(uh, npf*nfe, ncu);            
        // print2darray(q, npe, ncu);            
        // error("here");

        PGEMNMStridedBached(common.cublasHandle, npe, ncu, npe, one, &Cy[npe*npe*e1], npe, u, npe, zero, q, npe, ns, backend);
        PGEMNMStridedBached(common.cublasHandle, npe, ncu, npf*nfe, minusone, &Ey[npe*npf*nfe*e1], npe, uh, ndf, one, q, npe, ns, backend);
        if (common.wave==1) {
          ArrayExtract(s, sol.sdg, npe, nc, ne, 0, npe, 2*ncu, 3*ncu, e1, e2);          
          ArrayAXPBY(q, q, s, scalar, scalar, npe*ncu*ns); // scalar*(q + s)
        }
        ArrayInsert(udg, q, npe, nc, ne, 0, npe, 2*ncu, 3*ncu, e1, e2);
        //print3darray(q, npe, common.ncu, ns);         
      }
      else if (nd==3) {
        dstype *Cy = &res.C[npe*npe*common.ne];
        dstype *Cz = &res.C[npe*npe*common.ne*2];
        dstype *Ey = &res.E[npe*npf*nfe*common.ne];
        dstype *Ez = &res.E[npe*npf*nfe*common.ne*2];

        PGEMNMStridedBached(common.cublasHandle, npe, ncu, npe, one, &Cx[npe*npe*e1], npe, u, npe, zero, q, npe, ns, backend);
        PGEMNMStridedBached(common.cublasHandle, npe, ncu, npf*nfe, minusone, &Ex[npe*npf*nfe*e1], npe, uh, ndf, one, q, npe, ns, backend);
        if (common.wave==1) {
          ArrayExtract(s, sol.sdg, npe, nc, ne, 0, npe, ncu, 2*ncu, e1, e2);          
          ArrayAXPBY(q, q, s, scalar, scalar, npe*ncu*ns); // scalar*(q + s)
        }
        ArrayInsert(udg, q, npe, nc, ne, 0, npe, ncu, 2*ncu, e1, e2);

        PGEMNMStridedBached(common.cublasHandle, npe, ncu, npe, one, &Cy[npe*npe*e1], npe, u, npe, zero, q, npe, ns, backend);
        PGEMNMStridedBached(common.cublasHandle, npe, ncu, npf*nfe, minusone, &Ey[npe*npf*nfe*e1], npe, uh, ndf, one, q, npe, ns, backend);
        if (common.wave==1) {
          ArrayExtract(s, sol.sdg, npe, nc, ne, 0, npe, 2*ncu, 3*ncu, e1, e2);          
          ArrayAXPBY(q, q, s, scalar, scalar, npe*ncu*ns); // scalar*(q + s)
        }
        ArrayInsert(udg, q, npe, nc, ne, 0, npe, 2*ncu, 3*ncu, e1, e2);

        PGEMNMStridedBached(common.cublasHandle, npe, ncu, npe, one, &Cz[npe*npe*e1], npe, u, npe, zero, q, npe, ns, backend);                        
        PGEMNMStridedBached(common.cublasHandle, npe, ncu, npf*nfe, minusone, &Ez[npe*npf*nfe*e1], npe, uh, ndf, one, q, npe, ns, backend);
        if (common.wave==1) {
          ArrayExtract(s, sol.sdg, npe, nc, ne, 0, npe, 3*ncu, 4*ncu, e1, e2);          
          ArrayAXPBY(q, q, s, scalar, scalar, npe*ncu*ns); // scalar*(q + s)
        }
        ArrayInsert(udg, q, npe, nc, ne, 0, npe, 3*ncu, 4*ncu, e1, e2);
      }                
    }    
    
    if (common.debugMode==1) {
      string filename;
      if (common.mpiProcs==1)
        filename = common.fileout;
      else
        filename = common.fileout + NumberToString(common.mpiRank);      

      writearray2file(filename + "qEquation_uhat.bin", uhat, ncu*npf*common.nf, backend);  
      writearray2file(filename + "qEquation_udg.bin", udg, npe*nc*common.ne, backend);  
    }       
}

#endif

