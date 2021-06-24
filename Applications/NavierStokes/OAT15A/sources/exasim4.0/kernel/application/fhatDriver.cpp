#ifndef __FHATDRIVER
#define __FHATDRIVER

void FhatDriver(dstype *fg, dstype *xg, dstype *uh, dstype *ug1, dstype *ug2, dstype * og1, 
     dstype * og2, dstype *nl, meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, 
     tempstruct &tmp, commonstruct &common, Int ngf, Int f1, Int f2, Int backend)
{
    //Int nd = common.nd;     // spatial dimension       
    //Int ncu = common.ncu;// number of compoments of (u)    
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int ncq = common.ncq;// number of compoments of (q)
    Int ncp = common.ncp;// number of compoments of (p)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension            
    Int numPoints = ngf*(f2-f1);
    Int M = numPoints * ncu;
    Int N = numPoints*ncu*nd;
    Int ntau = common.ntau;
    dstype time = common.time; 
    
    // left flux
    FluxDriver(fg, xg, ug1, og1, mesh, master, app, sol, tmp, common,
                      ngf, f1, f2, backend);
    
    // right flux
    FluxDriver(&fg[N], xg, ug2, og2, mesh, master, app, sol, tmp, common,
                      ngf, f1, f2, backend);        
    
    // Part 1: fh = fg dot nl
    AverageFlux(fg, N, backend);    
    AverageFluxDotNormal(fg, nl, N, M, numPoints, nd, backend);            

    // Part 2: Contribution due to tau*(U-UH)    
    if (common.tau0<0) { 
#ifdef HAVE_ONETHREAD            
        if (backend==0) {
            opuStab(fg, xg, uh, ug1, ug2, og1, og2, nl, app.physicsparam, time, 
                             numPoints, nc, ncu, nd, ncx, nco, common.appname);
        }
#endif              
#ifdef HAVE_OPENMP                        
        if (backend==1) {
            cpuStab(fg, xg, uh, ug1, ug2, og1, og2, nl, app.physicsparam, time, 
                             numPoints, nc, ncu, nd, ncx, nco, common.appname);
        }
#endif            
#ifdef HAVE_CUDA                             
        if (backend==2) {
            gpuStab(fg, xg, uh, ug1, ug2, og1, og2, nl, app.physicsparam, time, 
                             numPoints, nc, ncu, nd, ncx, nco, common.appname);
        }
#endif                            
    }
    else if (ntau==0) {          
        AddStabilization1(fg, ug1, ug2, app.tau, M, backend);
    }
    else if (ntau==1) { // constant scalar  
//         if (backend==0) {            
//             ArrayCopy(&fg[N], &fg[0], N, backend);
//             print2darray(&fg[0], ngf, 4);
//             print2darray(&fg[N], ngf, 4);
//             opuStab(&fg[0], xg, uh, ug1, ug2, og1, og2, nl, app.physicsparam, time, 
//                              numPoints, nc, ncu, nd, ncx, nco, common.appname);
//         }            
         AddStabilization1(fg, ug1, ug2, app.tau, M, backend);
//         print2darray(&fg[0], ngf, 4);
//         print2darray(&fg[N], ngf, 4);        
    }
    else if (ntau==ncu) { // constant diagonal tensor
        AddStabilization2(fg, ug1, ug2, app.tau, M, numPoints, backend);
    }
    else if (ntau==ncu*ncu) { // constant full tensor      
        AddStabilization3(fg, ug1, ug2, app.tau, M, numPoints, ncu, backend);
    }
    else {
        printf("Stabilization option is not implemented");
        exit(-1);
    }            
}

void FhatDriver2(dstype *fg, dstype *xg, dstype *uh, dstype *ug1, dstype *ug2, dstype * og1, 
     dstype * og2, dstype *nl, meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, 
     tempstruct &tmp, commonstruct &common, Int ngf, Int f1, Int f2, Int backend)
{
    Int nd = common.nd;     // spatial dimension       
    Int ncu = common.ncu;// number of compoments of (u)
    Int numPoints = ngf*(f2-f1);
    Int M = numPoints * ncu;
    Int N = numPoints*ncu*nd;
    Int ntau = common.ntau;
    
    // Part 1: Contribution due to fh*nl    
    FluxDriver(&fg[N], xg, ug1, og1, mesh, master, app, sol, tmp, common,
                      ngf, f1, f2, backend);
    
    AverageFluxDotNormal(fg, nl, N, M, numPoints, nd, backend);

    // Part 2: Contribution due to tau*(U-UH)
    if (ntau==0) { 
        AddStabilization1(fg, ug1, uh, app.tau, M, backend);
    }
    else if (ntau==1) { // constant scalar  
        AddStabilization1(fg, ug1, uh, app.tau, M, backend);
    }
    else if (ntau==ncu) { // constant diagonal tensor
        AddStabilization2(fg, ug1, uh, app.tau, M, numPoints, backend);
    }
    else if (ntau==ncu*ncu) { // constant full tensor      
        AddStabilization3(fg, ug1, uh, app.tau, M, numPoints, ncu, backend);
    }
    else {
        printf("Stabilization option is not implemented");
        exit(-1);
    }
}

// void FhatDriver(dstype *fg, dstype *xg, dstype *uh, dstype *ug1, dstype *ug2, dstype * og1, 
//      dstype * og2, dstype *nl, meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, 
//      tempstruct &tmp, commonstruct &common, Int ngf, Int f1, Int f2, Int backend)
// {
//     Int nd = common.nd;     // spatial dimension       
//     Int ncu = common.ncu;// number of compoments of (u)
//     Int numPoints = ngf*(f2-f1);
//     Int M = numPoints * ncu;
//     Int N = numPoints*ncu*nd;
//     Int ntau = common.ntau;
//     
//     FluxDriver(fg, xg, ug1, og1, mesh, master, app, sol, tmp, common,
//                       ngf, f1, f2, backend);
//     FluxDriver(&fg[N], xg, ug2, og2, mesh, master, app, sol, tmp, common,
//                       ngf, f1, f2, backend);        
//     
// #ifdef HAVE_ONETHREAD           
//     if (backend==0) {
//         opuAverageFlux(fg, N);    
//         opuAverageFluxDotNormal(fg, nl, N, M, numPoints, nd);            
// 
//         // Part 2: Contribution due to tau*(U-UH)
//         if (ntau==0) { 
//             opuAddStabilization1(fg, ug1, ug2, app.tau, M);
//         }
//         else if (ntau==1) { // constant scalar  
//             opuAddStabilization1(fg, ug1, ug2, app.tau, M);
//         }
//         else if (ntau==ncu) { // constant diagonal tensor
//             opuAddStabilization2(fg, ug1, ug2, app.tau, M, numPoints);
//         }
//         else if (ntau==ncu*ncu) { // constant full tensor      
//             opuAddStabilization3(fg, ug1, ug2, app.tau, M, numPoints, ncu);
//         }
//         else {
//             printf("Stabilization option is not implemented");
//             exit(-1);
//         }            
//     }
// #endif              
// #ifdef HAVE_OPENMP            
//     if (backend==1) {
//         cpuAverageFlux(fg, N);                    
//         cpuAverageFluxDotNormal(fg, nl, N, M, numPoints, nd);            
//         
//         // Part 2: Contribution due to tau*(U-UH)
//         if (ntau==0) { 
//             cpuAddStabilization1(fg, ug1, ug2, app.tau, M);
//         }
//         else if (ntau==1) { // constant scalar  
//             cpuAddStabilization1(fg, ug1, ug2, app.tau, M);
//         }
//         else if (ntau==ncu) { // constant diagonal tensor
//             cpuAddStabilization2(fg, ug1, ug2, app.tau, M, numPoints);
//         }
//         else if (ntau==ncu*ncu) { // constant full tensor      
//             cpuAddStabilization3(fg, ug1, ug2, app.tau, M, numPoints, ncu);
//         }
//         else {
//             printf("Stabilization option is not implemented");
//             exit(-1);
//         }            
//     }
// #endif            
// #ifdef HAVE_CUDA                 
//     if (backend==2) {
//         gpuAverageFlux(fg, N);                            
//         gpuAverageFluxDotNormal(fg, nl, N, M, numPoints, nd);            
//         
//         // Part 2: Contribution due to tau*(U-UH)
//         if (ntau==0) { 
//             gpuAddStabilization1(fg, ug1, ug2, app.tau, M);
//         }
//         else if (ntau==1) { // constant scalar  
//             gpuAddStabilization1(fg, ug1, ug2, app.tau, M);
//         }
//         else if (ntau==ncu) { // constant diagonal tensor
//             gpuAddStabilization2(fg, ug1, ug2, app.tau, M, numPoints);
//         }
//         else if (ntau==ncu*ncu) { // constant full tensor      
//             gpuAddStabilization3(fg, ug1, ug2, app.tau, M, numPoints, ncu);
//         }
//         else {
//             printf("Stabilization option is not implemented");
//             exit(-1);
//         }            
//     }
// #endif        
// }
// 
// void FhatDriver2(dstype *fg, dstype *xg, dstype *uh, dstype *ug1, dstype *ug2, dstype * og1, 
//      dstype * og2, dstype *nl, meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, 
//      tempstruct &tmp, commonstruct &common, Int ngf, Int f1, Int f2, Int backend)
// {
//     Int nd = common.nd;     // spatial dimension       
//     Int ncu = common.ncu;// number of compoments of (u)
//     Int numPoints = ngf*(f2-f1);
//     Int M = numPoints * ncu;
//     Int N = numPoints*ncu*nd;
//     Int ntau = common.ntau;
//     
//     // Part 1: Contribution due to fh*nl    
//     FluxDriver(&fg[N], xg, ug1, og1, mesh, master, app, sol, tmp, common,
//                       ngf, f1, f2, backend);
//     
// #ifdef HAVE_ONETHREAD          
//     if (backend==0) {
//         opuAverageFluxDotNormal(fg, nl, N, M, numPoints, nd);
// 
//         // Part 2: Contribution due to tau*(U-UH)
//         if (ntau==0) { 
//             opuAddStabilization1(fg, ug1, uh, app.tau, M);
//         }
//         else if (ntau==1) { // constant scalar  
//             opuAddStabilization1(fg, ug1, uh, app.tau, M);
//         }
//         else if (ntau==ncu) { // constant diagonal tensor
//             opuAddStabilization2(fg, ug1, uh, app.tau, M, numPoints);
//         }
//         else if (ntau==ncu*ncu) { // constant full tensor      
//             opuAddStabilization3(fg, ug1, uh, app.tau, M, numPoints, ncu);
//         }
//         else {
//             printf("Stabilization option is not implemented");
//             exit(-1);
//         }
//     }
// #endif              
// #ifdef HAVE_OPENMP            
//     if (backend==1) {
//         cpuAverageFluxDotNormal(fg, nl, N, M, numPoints, nd);
// 
//         // Part 2: Contribution due to tau*(U-UH)
//         if (ntau==0) { 
//             cpuAddStabilization1(fg, ug1, uh, app.tau, M);
//         }
//         else if (ntau==1) { // constant scalar  
//             cpuAddStabilization1(fg, ug1, uh, app.tau, M);
//         }
//         else if (ntau==ncu) { // constant diagonal tensor
//             cpuAddStabilization2(fg, ug1, uh, app.tau, M, numPoints);
//         }
//         else if (ntau==ncu*ncu) { // constant full tensor      
//             cpuAddStabilization3(fg, ug1, uh, app.tau, M, numPoints, ncu);
//         }
//         else {
//             printf("Stabilization option is not implemented");
//             exit(-1);
//         }
//     }
// #endif            
// #ifdef HAVE_CUDA                 
//     if (backend==2) {
//         gpuAverageFluxDotNormal(fg, nl, N, M, numPoints, nd);
// 
//         // Part 2: Contribution due to tau*(U-UH)
//         if (ntau==0) { 
//             gpuAddStabilization1(fg, ug1, uh, app.tau, M);
//         }
//         else if (ntau==1) { // constant scalar  
//             gpuAddStabilization1(fg, ug1, uh, app.tau, M);
//         }
//         else if (ntau==ncu) { // constant diagonal tensor
//             gpuAddStabilization2(fg, ug1, uh, app.tau, M, numPoints);
//         }
//         else if (ntau==ncu*ncu) { // constant full tensor      
//             gpuAddStabilization3(fg, ug1, uh, app.tau, M, numPoints, ncu);
//         }
//         else {
//             printf("Stabilization option is not implemented");
//             exit(-1);
//         }            
//     }
// #endif        
// }
// 

#endif
