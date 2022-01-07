int IsElemCurved(dstype *jac, int ns, int nge)
{
    int curvedmesh = 0;
    for (int i=0; i<ns; i++) {
#ifdef HAVE_OPENMP                
        dstype minjac = cpuArrayMin(&jac[i*nge], nge);
        dstype maxjac = cpuArrayMax(&jac[i*nge], nge);  
#else                
        dstype minjac = opuArrayMin(&jac[i*nge], nge);
        dstype maxjac = opuArrayMax(&jac[i*nge], nge);          
#endif        
        if ((maxjac-minjac)>1e-10) {
            curvedmesh = 1;
            break;
        }
    }
    return curvedmesh;
}

int IsMeshCurved(solstruct &sol, appstruct &app, masterstruct &master, meshstruct &mesh, tempstruct &tmp)
{            
    dstype *xn, *Xx, *jac, *Jg;  
    
    Int ncx = app.ndims[11];// number of compoments of (xdg)        
    Int nd = master.ndims[0];     // spatial dimension    
    Int npe = master.ndims[5]; // number of nodes on master element
    Int nge = master.ndims[7]; // number of gauss points on master element    
    //Int ne = mesh.ndims[1]; // number of elements in this subdomain 
    Int nbe = mesh.ndims[5]; // number of blocks for elements            
    
    // compute volume integrals to form Rq
    int curvedmesh = 0;
    for (Int j=0; j<nbe; j++) // for each block of elements
    {
        Int e1 = mesh.eblks[3*j]-1;
        Int e2 = mesh.eblks[3*j+1];
        Int ns = e2-e1;        
        Int nga = nge*ns;
        
        xn = &tmp.tempn[0];
        Xx = &tmp.tempg[0];
        jac = &tmp.tempg[nga*(nd*nd)];
        Jg = &tmp.tempg[nga*(nd*nd+1)];
                
#ifdef HAVE_OPENMP        
        cpuGetElemNodes(xn, sol.xdg, npe, ncx, 0, ncx, e1, e2);                                
        for (Int k=0; k<nd; k++)
            cpuNode2Gauss(&Jg[k*nga*nd], xn, &master.shapegt[(k+1)*nge*npe], nge, npe, ns*nd);                 
        cpuElemGeom(Xx, jac, Jg, ns, nge, nd);
#else        
        opuGetElemNodes(xn, sol.xdg, npe, ncx, 0, ncx, e1, e2);                                
        for (Int k=0; k<nd; k++)
            cpuNode2Gauss(&Jg[k*nga*nd], xn, &master.shapegt[(k+1)*nge*npe], nge, npe, ns*nd);                 
        opuElemGeom(Xx, jac, Jg, ns, nge, nd);        
#endif        
        if (IsElemCurved(jac, ns, nge)) {
            curvedmesh = 1;
            break;
        }        
    }                
    
    return curvedmesh;
}

