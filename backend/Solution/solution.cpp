#ifndef __SOLUTION
#define __SOLUTION

#include "solution.h"
#include "previoussolutions.cpp"
#include "updatesolution.cpp"
#include "updatesource.cpp"
#include "timestepcoeff.cpp"
#include "avsolution.cpp"

#ifdef TIMESTEP  
#include <sys/time.h>
#endif

void CSolution::SteadyProblem(ofstream &out, Int backend) 
{   
    INIT_TIMING;        
#ifdef TIMING    
    for (int i=0; i<100; i++)
        disc.common.timing[i] = 0.0; 
#endif
    // obtain odg from the solutions of the other PDE models
    if (disc.common.nomodels>1) {
        Int nco = disc.common.nco;
        Int npe = disc.common.npe;
        Int ne = disc.common.ne;            
        for (Int n=0; n<nco; n++) {            
            Int m = disc.common.vindx[n];     // model index
            Int k = disc.common.vindx[nco+n]; // solution index
            // extract the kth component of udg from PDE model m and store it in Ru
            ArrayExtract(disc.res.Ru, disc.sol.udgarray[m], npe, disc.common.ncarray[m], ne, 0, npe, k, k+1, 0, ne);         
            // insert Ru into odg
            ArrayInsert(disc.sol.odg, disc.res.Ru, npe, nco, ne, 0, npe, n, n+1, 0, ne);          
        }
    }
    
    // calculate AV field
    if (disc.common.ncAV>0 && disc.common.frozenAVflag > 0) {
        // START_TIMING;
        // compute AV field at sys.u and store it in sol.odg
        //disc.evalAVfield(disc.sol.odg, solv.sys.u, backend);
        disc.evalAVfield(disc.sol.odg, backend);

        for (Int iav = 0; iav<disc.common.AVsmoothingIter; iav++){
            // printf("Solution AV smoothing iter: %i\n", iav);
            disc.DG2CG2(disc.sol.odg, disc.sol.odg, solv.sys.x, disc.common.nco, disc.common.nco, disc.common.ncAV, backend);

#ifdef  HAVE_MPI    
            Int bsz = disc.common.npe*disc.common.ncAV;
            Int nudg = disc.common.npe*disc.common.nco;
            Int n;

            //for (n=0; n<disc.common.nelemsend; n++)
            //    ArrayCopy(&disc.tmp.buffsend[bsz*n], &disc.sol.odg[nudg*disc.common.elemsend[n]], bsz, backend);
            GetArrayAtIndex(disc.tmp.buffsend, disc.sol.odg, disc.mesh.elemsendodg, bsz*disc.common.nelemsend);

// #ifdef HAVE_CUDA
//             cudaDeviceSynchronize();
// #endif

#ifdef HAVE_HIP
            hipDeviceSynchronize();
#endif
            
            Int neighbor, nsend, psend = 0, request_counter = 0;
            for (n=0; n<disc.common.nnbsd; n++) {
                neighbor = disc.common.nbsd[n];
                nsend = disc.common.elemsendpts[n]*bsz;
                if (nsend>0) {
                    MPI_Isend(&disc.tmp.buffsend[psend], nsend, MPI_DOUBLE, neighbor, 0,
                        MPI_COMM_WORLD, &disc.common.requests[request_counter]);
                    psend += nsend;
                    request_counter += 1;
                }
            }

            Int nrecv, precv = 0;
            for (n=0; n<disc.common.nnbsd; n++) {
                neighbor = disc.common.nbsd[n];
                nrecv = disc.common.elemrecvpts[n]*bsz;
                if (nrecv>0) {
                    MPI_Irecv(&disc.tmp.buffrecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                        MPI_COMM_WORLD, &disc.common.requests[request_counter]);
                    precv += nrecv;
                    request_counter += 1;
                }
            }

            MPI_Waitall(request_counter, disc.common.requests, disc.common.statuses);
            //for (n=0; n<disc.common.nelemrecv; n++)
            //   ArrayCopy(&disc.sol.odg[nudg*disc.common.elemrecv[n]], &disc.tmp.buffrecv[bsz*n], bsz, backend);
            PutArrayAtIndex(disc.sol.odg, disc.tmp.buffrecv, disc.mesh.elemrecvodg, bsz*disc.common.nelemrecv);
#endif
    //    END_TIMING_DISC(98);    
        }

    }
    if (disc.common.nco>0) {
        for (Int j=0; j<disc.common.nbe; j++) {
            Int e1 = disc.common.eblks[3*j]-1;
            Int e2 = disc.common.eblks[3*j+1];                
            GetElemNodes(disc.tmp.tempn, disc.sol.odg, disc.common.npe, disc.common.nco, 
                    0, disc.common.nco, e1, e2);        
            Node2Gauss(disc.common.cublasHandle, &disc.sol.odgg[disc.common.nge*disc.common.nco*e1], 
              disc.tmp.tempn, disc.master.shapegt, disc.common.nge, disc.common.npe, (e2-e1)*disc.common.nco, backend);        
        }         
        for (Int j=0; j<disc.common.nbf; j++) {
            Int f1 = disc.common.fblks[3*j]-1;
            Int f2 = disc.common.fblks[3*j+1];            
            
            GetFaceNodes(disc.tmp.tempn, disc.sol.odg, disc.mesh.facecon, disc.common.npf, disc.common.nco, 
                    disc.common.npe, disc.common.nco, f1, f2, 1);          
            Node2Gauss(disc.common.cublasHandle, &disc.sol.og1[disc.common.ngf*disc.common.nco*f1], 
              disc.tmp.tempn, disc.master.shapfgt, disc.common.ngf, disc.common.npf, (f2-f1)*disc.common.nco, backend);               
            
            GetFaceNodes(disc.tmp.tempn, disc.sol.odg, disc.mesh.facecon, disc.common.npf, disc.common.nco, 
                    disc.common.npe, disc.common.nco, f1, f2, 2);          
            Node2Gauss(disc.common.cublasHandle, &disc.sol.og2[disc.common.ngf*disc.common.nco*f1], 
              disc.tmp.tempn, disc.master.shapfgt, disc.common.ngf, disc.common.npf, (f2-f1)*disc.common.nco, backend);               
        }        
    }
    
    if (disc.common.ncs>0) {
        for (Int j=0; j<disc.common.nbe; j++) {
            Int e1 = disc.common.eblks[3*j]-1;
            Int e2 = disc.common.eblks[3*j+1];                
            GetElemNodes(disc.tmp.tempn, disc.sol.sdg, disc.common.npe, disc.common.ncs, 0, disc.common.ncs, e1, e2);        
            Node2Gauss(disc.common.cublasHandle, &disc.sol.sdgg[disc.common.nge*disc.common.ncs*e1], 
              disc.tmp.tempn, disc.master.shapegt, disc.common.nge, disc.common.npe, (e2-e1)*disc.common.ncs, backend);        
        } 
    }
    
    // use PTC to solve steady problem
    if (disc.common.spatialScheme==0) {
      solv.PseudoTransientContinuation(disc, prec, out, backend);           
    }
    else if (disc.common.spatialScheme==1) {      
      solv.NewtonSolver(disc, prec, out, disc.common.ndofuhat, disc.common.spatialScheme, backend);           
    }
    else
      error("Spatial discretization scheme is not implemented");
        
#ifdef TIMING         
    if (disc.common.mpiRank==0) {
        printf("\nComputing initial guess time: %g miliseconds\n", disc.common.timing[99]);   
        printf("Computing AV fields time: %g miliseconds\n", disc.common.timing[98]);   
        printf("Nonlinear solver time: %g miliseconds\n", disc.common.timing[97]);                
    }
#endif    
}

void CSolution::InitSolution(Int backend) 
{    
//     // compute the geometry quantities
//     disc.compGeometry(backend);
//     
//     //printArray2D(&disc.sol.elemg[0],disc.common.nge,10,backend);
//     
//     // compute the inverse of the mass matrix
//     disc.compMassInverse(backend);
            
    // compute q
    // if ((disc.common.ncq>0) & (disc.common.wave==0))
    //     disc.evalQSer(backend);
            
    // // set pointer depending on the matrix type
    // if (disc.common.precMatrixType==0)
    //     prec.precond.Cmat = &prec.precond.C[0];
    // else //if (disc.common.precMatrixType==2)
    //     prec.precond.Cmat = &disc.res.Minv[0];    
        
    if (disc.common.spatialScheme==0) {
        ArrayExtract(solv.sys.u, disc.sol.udg, disc.common.npe, disc.common.nc, disc.common.ne1, 
              0, disc.common.npe, 0, disc.common.ncu, 0, disc.common.ne1);                                                  
    }
    else if (disc.common.spatialScheme==1) {      
        ArrayCopy(solv.sys.u, disc.sol.uh, disc.common.ndofuhat);
    }
    else
        error("Spatial discretization scheme is not implemented");
            
    // save solutions into binary files
    this->SaveNodesOnBoundary(backend);     
    
    if (disc.common.tdep==1) { // DIRK schemes
        //DIRK coefficients 
        disc.common.temporalScheme = 0; 
        TimestepCoefficents(disc.common); 
                
        if (disc.common.mpiRank==0)
            cout<<"Compute solution average = "<<disc.common.compudgavg<<endl;
        
        if (disc.common.compudgavg == 1) {
            string filename = disc.common.fileout + "avg_np" + NumberToString(disc.common.mpiRank) + ".bin";
            disc.common.readudgavg = fileexists(filename);
            if (disc.common.mpiRank==0)
                cout<<"File exist = "<<disc.common.readudgavg<<endl;
            if (disc.common.readudgavg == 0)
                ArraySetValue(disc.sol.udgavg, zero, disc.common.ndofudg1+1);
            else 
                readarrayfromfile(filename, &disc.sol.udgavg, disc.common.ndofudg1+1, backend);   
        }        
    }    
}

void CSolution::DIRK(ofstream &out, Int backend)
{    
    INIT_TIMING;        
    
    // initial time
    dstype time = disc.common.time;           
    
    //DIRK coefficients 
    disc.common.temporalScheme = 0; 
    TimestepCoefficents(disc.common); 
                
#ifdef TIMESTEP                  
    struct timeval tv1, tv2;
#endif                
    
    // time stepping with DIRK schemes
    for (Int istep=0; istep<disc.common.tsteps; istep++)            
    {            
        // current timestep        
        disc.common.currentstep = istep;
        
        // store previous solutions to calculate the source term        
        PreviousSolutions(disc.sol, solv.sys, disc.common, backend);

#ifdef TIMESTEP              
        gettimeofday(&tv1, NULL); 
#endif
                    
        // compute the solution at the next step
        for (Int j=0; j<disc.common.tstages; j++) {            
            // current timestage
            disc.common.currentstage = j;
        
            // current time
            disc.common.time = time + disc.common.dt[istep]*disc.common.DIRKcoeff_t[j];

            if (disc.common.mpiRank==0)
                printf("\nTimestep :  %d,  Timestage :  %d,   Time : %g\n",istep+1,j+1,disc.common.time);            
        
#ifdef TIMING    
            disc.common.timing[100] = 0.0; 
            disc.common.timing[101] = 0.0; 
#endif
        
            START_TIMING;

            // update source term             
            UpdateSource(disc.sol, solv.sys, disc.app, disc.res, disc.common, backend);
            END_TIMING_DISC(100);    

            // solve the problem 
            this->SteadyProblem(out, backend);                             

            START_TIMING;
            // update solution 
            UpdateSolution(disc.sol, solv.sys, disc.common, backend);                     
            END_TIMING_DISC(101);

#ifdef TIMING         
            if (disc.common.mpiRank==0) {
                printf("Updating source term time: %g miliseconds\n", disc.common.timing[100]);   
                printf("Updating solution time: %g miliseconds\n\n", disc.common.timing[101]);                           
            }
#endif            
        }
        
        //compute time-average solution
        if (disc.common.compudgavg == 1) {
            ArrayAXPBY(disc.sol.udgavg, disc.sol.udgavg, disc.sol.udg, one, one, disc.common.ndofudg1);            
            ArrayAddScalar(&disc.sol.udgavg[disc.common.ndofudg1], one, 1);
        }

        // save solutions into binary files
        //SaveSolutions(disc.sol, solv.sys, disc.common, backend);            
        this->SaveSolutions(backend); 
        this->SaveSolutionsOnBoundary(backend); 
        if (disc.common.nce>0)
            this->SaveOutputCG(backend);    
        
#ifdef TIMESTEP                          
        gettimeofday(&tv2, NULL);            
        if (disc.common.mpiRank==0)
            printf("\nExecution time (in millisec) for timestep %d:  %g\n", istep+1,
                (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
                (double)(tv2.tv_sec -tv1.tv_sec )*1000);
#endif                    
        // update time
        time = time + disc.common.dt[istep];                    
    }           
}

void CSolution::SteadyProblem_PTC(ofstream &out, Int backend) {

    // initial time
    double time = disc.common.time;           
    double monitor_diff, monitor_scale, delta_monitor;
    int N = disc.common.ndofuhat;
    int NLiters = disc.common.nonlinearSolverMaxIter;
    double nrmr = 0;
    int conv_flag = 0;

    // time stepping with DIRK schemes
    for (Int istep=0; istep<disc.common.tsteps-1; istep++)            
    {            
        disc.common.nonlinearSolverMaxIter = 1;

        // current timestep        
        disc.common.currentstep = istep;

        // store previous solutions to calculate the source term        
        PreviousSolutions(disc.sol, solv.sys, disc.common, backend);
                            
        // compute the solution at the next step
        for (int j=0; j<disc.common.tstages; j++) {     
            
            if (disc.common.mpiRank==0)
                printf("\nTimestep :  %d,  Timestage :  %d,   Time : %g\n",istep+1,j+1,time + disc.common.dt[istep]*disc.common.DIRKcoeff_t[j]);                                
                            
            // current timestage
            disc.common.currentstage = j;

            // current time
            disc.common.time = time + disc.common.dt[istep]*disc.common.DIRKcoeff_t[j];

            // update source term             
            UpdateSource(disc.sol, solv.sys, disc.app, disc.res, disc.common, backend);
            
            // solve the problem 
            this->SteadyProblem(out, backend);                             

            // update solution 
            UpdateSolution(disc.sol, solv.sys, disc.app, disc.res, disc.tmp, disc.common, backend);
            
            // TODO: input wprev
            disc.evalMonitor(disc.tmp.tempn,  disc.sol.udg, disc.sol.wdg, disc.common.nc, backend);
            disc.evalMonitor(disc.tmp.tempg,  solv.sys.udgprev, solv.sys.wprev, disc.common.ncu, backend);
            
            ArrayAXPBY(disc.tmp.tempn, disc.tmp.tempn, disc.tmp.tempg, 1.0, -1.0, disc.common.npe*disc.common.ncm*disc.common.ne);            
            
            monitor_diff  = PNORM(disc.common.cublasHandle,  disc.common.npe*disc.common.ncm*disc.common.ne,disc.tmp.tempn, backend);
            monitor_scale = PNORM(disc.common.cublasHandle,  disc.common.npe*disc.common.ncm*disc.common.ne, disc.tmp.tempg, backend);

            delta_monitor = monitor_diff / monitor_scale;
            std::cout << "delta_monitor: " << delta_monitor << std::endl;

            if ((delta_monitor > 1.0 || solv.sys.alpha < 0.1))
            {
                std::cout << "Linesearch failed or large change in solution: reducing timestep" << std::endl;
                
                // Revert time step
                disc.common.time = disc.common.time - disc.common.dt[istep]*disc.common.DIRKcoeff_t[j];

                // Copy udg old to udg
                ArrayCopy(disc.sol.udg, solv.sys.udgprev, disc.common.npe*disc.common.ncu*disc.common.ne2);

                // Compute a new UH
                GetFaceNodes(disc.sol.uh, disc.sol.udg, disc.mesh.f2e, disc.mesh.perm, disc.common.npf, disc.common.ncu, disc.common.npe, disc.common.nc, disc.common.nf);

                // Recompute gradient from udg old
                hdgGetQ(disc.sol.udg, disc.sol.uh, disc.sol, disc.res, disc.mesh, disc.tmp, disc.common, backend);

                // decrease timestep by 10
                disc.common.dt[istep+1] = disc.common.dt[istep] / 10;
                std::cout << "next time step: " << disc.common.dt[istep+1] << std::endl;
                if (disc.common.dt[istep+1] < 1e-8){
                    std::cout << "WARNING: PTC stalled" << std::endl;
                    istep = disc.common.tsteps+1;
                }
            }
            else if (delta_monitor < 0.1)
            {
                if (disc.common.linearSolverIter < disc.common.linearSolverMaxIter){
                    // increase timestep by 2
                    disc.common.dt[istep+1] = disc.common.dt[istep]*2;
                }
                else{
                    std::cout << "Too many GMRES iterations, not increasing timestep" << std::endl;
                    disc.common.dt[istep+1] = disc.common.dt[istep]*1;
                }
                

                if (delta_monitor < 1e-4 && disc.common.dt[istep] > 1e-3) {
                    if (disc.common.runmode == 10) {
                        std::cout << "Evaluate steady residual..." << std::endl;
                        disc.common.tdep=0;
    
                        if (disc.common.ncq > 0) hdgGetQ(disc.sol.udg, disc.sol.uh, disc.sol, disc.res, disc.mesh, disc.tmp, disc.common, backend);          
            
                        // compute the residual vector R = [Ru; Rh]
                        disc.hdgAssembleResidual(solv.sys.b, backend);
                                
                        nrmr = PNORM(disc.common.cublasHandle, N, solv.sys.b, backend);       
                        nrmr += PNORM(disc.common.cublasHandle, disc.common.npe*disc.common.ncu*disc.common.ne1, disc.res.Ru, backend); 
                        std::cout << " Steady residual = " << nrmr << std::endl;
    
                        // SaveSolutions(backend); 
                        if (nrmr < disc.common.nonlinearSolverTol) {
                            conv_flag = 1;
                            istep = disc.common.tsteps+10;
                            this->SaveSolutions(backend); 
                            this->SaveSolutionsOnBoundary(backend); 
                        }
                        // istep = disc.common.tsteps+1;
                        disc.common.tdep=1;
                    }
                    if (disc.common.runmode == 11) { // Compute steady solve
                        std::cout << "Steady solve..." << std::endl;
                        disc.common.tdep=0;
                        disc.common.nonlinearSolverMaxIter = NLiters;
                        this->SolveProblem(out, backend);
                        istep = disc.common.tsteps+10;
                    }
                    
                }    
            }
        }
        // update time
        time = time + disc.common.dt[istep];                    
    }   
    if (conv_flag == 0) {                
        std::cout << "Warning: PTC reached max iterations without converging." << std::endl;
        // Save steady solution anyways
        disc.common.tdep=0;
        this->SaveSolutions(backend); 
        this->SaveSolutionsOnBoundary(backend); 
    }
}

void CSolution::SolveProblem(ofstream &out, Int backend) 
{          
    this->InitSolution(backend); 
        
    if (disc.common.tdep==1) {        
        // solve time-dependent problems using DIRK
        this->DIRK(out, backend);            
    }
    else {
        // solve steady-state problems
        this->SteadyProblem(out, backend);        
                
        // save solutions into binary files            
        this->SaveSolutions(backend);    
                
        this->SaveSolutionsOnBoundary(backend);         
        if (disc.common.nce>0)
            this->SaveOutputCG(backend);            
    }        
}

void CSolution::SaveSolutions(Int backend) 
{
   if (disc.common.tdep==1) { 
        if (((disc.common.currentstep+1) % disc.common.saveSolFreq) == 0)             
        {        
            string filename = disc.common.fileout + "_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";     
            if (disc.common.saveSolOpt==0)
                writearray2file(filename, solv.sys.u, disc.common.ndof1, backend);
            else
                writearray2file(filename, disc.sol.udg, disc.common.ndofudg1, backend);
                        
            if (disc.common.ncw>0) {
                string fn = disc.common.fileout + "_wdg_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";                    
                writearray2file(fn, solv.sys.wtmp, disc.common.ndofw1, backend);
            }                        
            
            if (disc.common.compudgavg == 1) {
                string fn1 = disc.common.fileout + "avg_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin"; 
                writearray2file(fn1, disc.sol.udgavg, disc.common.ndofudg1+1, backend);
            }

            if (disc.common.spatialScheme==1) {
                string fn2 = disc.common.fileout + "_uhat_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";                    
                writearray2file(fn2, disc.sol.uh, disc.common.ndofuhat, backend);        
            }
        }    
   }
   else {
        string filename = disc.common.fileout + "_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";                    
        if (disc.common.saveSolOpt==0)
            writearray2file(filename, solv.sys.u, disc.common.ndof1, backend);
        else
            writearray2file(filename, disc.sol.udg, disc.common.ndofudg1, backend);       
        
        if (disc.common.ncw>0) {
            string fn = disc.common.fileout + "_wdg_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";                    
            writearray2file(fn, disc.sol.wdg, disc.common.ndofw1, backend);     
        }                

        if (disc.common.spatialScheme==1) {
            string filename1 = disc.common.fileout + "_uhat_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";                    
            writearray2file(filename1, disc.sol.uh, disc.common.ndofuhat, backend);        
        }
   }    
}

void CSolution::ReadSolutions(Int backend) 
{
   if (disc.common.tdep==1) { 
        if (((disc.common.currentstep+1) % disc.common.saveSolFreq) == 0)             
        {        
            string filename = disc.common.fileout + "_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";     
            if (disc.common.saveSolOpt==0) {
                readarrayfromfile(filename, &solv.sys.u, disc.common.ndof1, backend);
                // insert u into udg
                ArrayInsert(disc.sol.udg, solv.sys.u, disc.common.npe, disc.common.nc, 
                 disc.common.ne, 0, disc.common.npe, 0, disc.common.ncu, 0, disc.common.ne1);  
            }
            else
                readarrayfromfile(filename, &disc.sol.udg, disc.common.ndofudg1, backend);        
            
            if (disc.common.ncw>0) {
                string fn = disc.common.fileout+"_wdg_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";                    
                readarrayfromfile(fn, &disc.sol.wdg, disc.common.ndofw1, backend);     
            }                      

            if (disc.common.spatialScheme==1) {
                string fn2 = disc.common.fileout + "_uhat_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";                    
                readarrayfromfile(fn2, &disc.sol.uh, disc.common.ndofuhat, backend);        
            }              
        }                                
   }
   else {
        string filename = disc.common.fileout + "_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";                    
        if (disc.common.saveSolOpt==0) {
            readarrayfromfile(filename, &solv.sys.u, disc.common.ndof1, backend);
            // insert u into udg
            ArrayInsert(disc.sol.udg, solv.sys.u, disc.common.npe, disc.common.nc, 
             disc.common.ne, 0, disc.common.npe, 0, disc.common.ncu, 0, disc.common.ne1);              
        }
        else
            readarrayfromfile(filename, &disc.sol.udg, disc.common.ndofudg1, backend);      
        
        if (disc.common.ncw>0) {
            string fn = disc.common.fileout + "_wdg_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";                    
            readarrayfromfile(fn, &disc.sol.wdg, disc.common.ndofw1, backend);     
        }                

        if (disc.common.spatialScheme==1) {
            string fn = disc.common.fileout + "_uhat_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";                    
            readarrayfromfile(fn, &disc.sol.uh, disc.common.ndofuhat, backend);        
        }                                    
   }    
}
 
void CSolution::SaveOutputDG(Int backend) 
{
   if (disc.common.tdep==1) { 
        if (((disc.common.currentstep+1) % disc.common.saveSolFreq) == 0)             
        {                    
            string filename1 = disc.common.fileout + "_outputDG_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";     
            disc.evalOutput(solv.sys.v, backend);                        
            writearray2file(filename1, solv.sys.v, disc.common.ndofedg1, backend);       
        }                                
   }
   else {
        string filename1 = disc.common.fileout + "_outputDG_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";                           
        disc.evalOutput(solv.sys.v, backend);
        writearray2file(filename1, solv.sys.v, disc.common.ndofedg1, backend);       
   }    
}

void CSolution::SaveOutputCG(Int backend) 
{
   if (disc.common.tdep==1) { 
        if (((disc.common.currentstep+1) % disc.common.saveSolFreq) == 0)             
        {                    
            string filename1 = disc.common.fileout + "_outputCG_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";     
            disc.evalOutput(solv.sys.v, backend);
            disc.DG2CG(solv.sys.v, solv.sys.v, solv.sys.x, disc.common.nce, 
                     disc.common.nce, disc.common.nce, backend);
            writearray2file(filename1, solv.sys.v, disc.common.ndofedg1, backend);                   
//             disc.DG2CG3(solv.sys.v, solv.sys.v, solv.sys.x, disc.common.nce, 
//                  disc.common.nce, disc.common.nce, backend);
//             writearray2file(filename1, solv.sys.v, disc.common.ndofucg, backend);               
        }                                
   }
   else {
        string filename1 = disc.common.fileout + "_outputCG_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";                            
        disc.evalOutput(solv.sys.v, backend);
        disc.DG2CG(solv.sys.v, solv.sys.v, solv.sys.x, disc.common.nce, 
                 disc.common.nce, disc.common.nce, backend);
        writearray2file(filename1, solv.sys.v, disc.common.ndofedg1, backend);               
//         disc.DG2CG3(solv.sys.v, solv.sys.v, solv.sys.x, disc.common.nce, 
//                  disc.common.nce, disc.common.nce, backend);        
//         writearray2file(filename1, solv.sys.v, disc.common.ndofucg, backend);               
   }    
}        

void CSolution::SaveSolutionsOnBoundary(Int backend) 
{   
    if ( disc.common.saveSolBouFreq>0 ) {
        if (((disc.common.currentstep+1) % disc.common.saveSolBouFreq) == 0)             
        {        
            for (Int j=0; j<disc.common.nbf; j++) {
                Int f1 = disc.common.fblks[3*j]-1;
                Int f2 = disc.common.fblks[3*j+1];    
                Int ib = disc.common.fblks[3*j+2];            
                if (ib == disc.common.ibs) {     
                    Int npf = disc.common.npf; // number of nodes on master face      
                    Int nf = f2-f1;
                    Int nn = npf*nf; 
                    Int nc = disc.common.nc; // number of compoments of (u, q, p)            
                    GetArrayAtIndex(disc.tmp.tempn, disc.sol.udg, &disc.mesh.findudg1[npf*nc*f1], nn*nc);
                    string filename = disc.common.fileout + "bou_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";     
                    writearray2file(filename, disc.tmp.tempn, nn*nc, backend);            
                }
            }                                               
        }                                
    }
}

void CSolution::SaveNodesOnBoundary(Int backend) 
{   
    if ( disc.common.saveSolBouFreq>0 ) {
        for (Int j=0; j<disc.common.nbf; j++) {
            Int f1 = disc.common.fblks[3*j]-1;
            Int f2 = disc.common.fblks[3*j+1];    
            Int ib = disc.common.fblks[3*j+2];            
            if (ib == disc.common.ibs) {     
                Int npf = disc.common.npf; // number of nodes on master face      
                Int nf = f2-f1;
                Int nn = npf*nf; 
                Int ncx = disc.common.ncx; // number of compoments of (u, q, p)                            
                GetArrayAtIndex(disc.tmp.tempn, disc.sol.xdg, &disc.mesh.findxdg1[npf*ncx*f1], nn*ncx);
                string filename = disc.common.fileout + "node_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";     
                writearray2file(filename, disc.tmp.tempn, nn*ncx, backend);            
            }
        }                                                                   
    }
}

#endif        
