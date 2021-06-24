#ifndef __SOLUTION
#define __SOLUTION

#include "solution.h"
#include "previoussolutions.cpp"
#include "savesolutions.cpp"
#include "updatesolution.cpp"
#include "updatesource.cpp"
#include "timestepcoeff.cpp"

#ifdef TIMESTEP  
#include <sys/time.h>
#endif

void CSolution::SteadyProblem(Int backend) 
{   
    INIT_TIMING;        
#ifdef TIMING    
    for (int i=0; i<100; i++)
        disc.common.timing[i] = 0.0; 
#endif
    
    // use minimal residual to compute the initial guess
//     if (disc.common.RBcurrentdim>0) {
//         START_TIMING;
//         solv.InitialGuess(disc, prec, backend);
//         END_TIMING_DISC(99);    
//     }
        
//     // calculate AV field
    if (disc.common.AVflag>0) {
        START_TIMING;
        // compute AV field at sys.u and store it in sol.odg
        //disc.evalAVfield(disc.sol.odg, solv.sys.u, backend);
        disc.evalAVfield(disc.sol.odg, backend);
        
        //writearray2file(disc.common.fileout + "udg" + "_np" + NumberToString(disc.common.mpiRank) + ".bin", disc.sol.odg, disc.common.npe*disc.common.nco*disc.common.ne1, backend);
        
        // smooth AV field
        disc.DG2CG2(disc.sol.odg, disc.sol.odg, solv.sys.x, disc.common.nco, 
               disc.common.nco, disc.common.AVflag, backend);

            
#ifdef  HAVE_MPI    
        Int bsz = disc.common.npe*disc.common.nco;
        Int nudg = disc.common.npe*disc.common.nco;
        Int n;

        /* copy some portion of u to buffsend */
        for (n=0; n<disc.common.nelemsend; n++)         
            ArrayCopy(&disc.tmp.buffsend[bsz*n], &disc.sol.odg[nudg*disc.common.elemsend[n]], bsz, backend);            

        /* non-blocking send */
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

        /* non-blocking receive */
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

        // non-blocking receive solutions on exterior and outer elements from neighbors
        /* wait until all send and receive operations are completely done */
        MPI_Waitall(request_counter, disc.common.requests, disc.common.statuses);

        /* copy buffrecv to udg */
        for (n=0; n<disc.common.nelemrecv; n++) 
            ArrayCopy(&disc.sol.odg[nudg*disc.common.elemrecv[n]], &disc.tmp.buffrecv[bsz*n], bsz, backend);        
#endif
             
        disc.DG2CG2(disc.sol.odg, disc.sol.odg, solv.sys.x, disc.common.nco, 
               disc.common.nco, disc.common.AVflag, backend);
                 
#ifdef  HAVE_MPI         
        for (n=0; n<disc.common.nelemsend; n++)         
            ArrayCopy(&disc.tmp.buffsend[bsz*n], &disc.sol.odg[nudg*disc.common.elemsend[n]], bsz, backend);            

        psend = 0;
        request_counter = 0;
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

        precv = 0;
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

        // non-blocking receive solutions on exterior and outer elements from neighbors
        MPI_Waitall(request_counter, disc.common.requests, disc.common.statuses);

        for (n=0; n<disc.common.nelemrecv; n++) 
            ArrayCopy(&disc.sol.odg[nudg*disc.common.elemrecv[n]], &disc.tmp.buffrecv[bsz*n], bsz, backend);        
#endif   
        END_TIMING_DISC(98);    
    }
    //error("here");
    
    if (disc.common.nco>0) {
        for (Int j=0; j<disc.common.nbe; j++) {
            Int e1 = disc.common.eblks[3*j]-1;
            Int e2 = disc.common.eblks[3*j+1];                
            GetElemNodes(disc.tmp.tempn, disc.sol.odg, disc.common.npe, disc.common.nco, 
                    0, disc.common.nco, e1, e2, backend);        
            Node2Gauss(disc.common.cublasHandle, &disc.sol.odgg[disc.common.nge*disc.common.nco*e1], 
              disc.tmp.tempn, disc.master.shapegt, disc.common.nge, disc.common.npe, (e2-e1)*disc.common.nco, backend);        
        }         
        for (Int j=0; j<disc.common.nbf; j++) {
            Int f1 = disc.common.fblks[3*j]-1;
            Int f2 = disc.common.fblks[3*j+1];            
            
            GetFaceNodes(disc.tmp.tempn, disc.sol.odg, disc.mesh.facecon, disc.common.npf, disc.common.nco, 
                    disc.common.npe, disc.common.nco, f1, f2, 1, backend);          
            Node2Gauss(disc.common.cublasHandle, &disc.sol.og1[disc.common.ngf*disc.common.nco*f1], 
              disc.tmp.tempn, disc.master.shapfgt, disc.common.ngf, disc.common.npf, (f2-f1)*disc.common.nco, backend);               
            
            GetFaceNodes(disc.tmp.tempn, disc.sol.odg, disc.mesh.facecon, disc.common.npf, disc.common.nco, 
                    disc.common.npe, disc.common.nco, f1, f2, 2, backend);          
            Node2Gauss(disc.common.cublasHandle, &disc.sol.og2[disc.common.ngf*disc.common.nco*f1], 
              disc.tmp.tempn, disc.master.shapfgt, disc.common.ngf, disc.common.npf, (f2-f1)*disc.common.nco, backend);               
        }        
    }
    
    if (disc.common.ncs>0) {
        for (Int j=0; j<disc.common.nbe; j++) {
            Int e1 = disc.common.eblks[3*j]-1;
            Int e2 = disc.common.eblks[3*j+1];                
            GetElemNodes(disc.tmp.tempn, disc.sol.sdg, disc.common.npe, disc.common.ncs, 0, disc.common.ncs, e1, e2, backend);        
            Node2Gauss(disc.common.cublasHandle, &disc.sol.sdgg[disc.common.nge*disc.common.ncs*e1], 
              disc.tmp.tempn, disc.master.shapegt, disc.common.nge, disc.common.npe, (e2-e1)*disc.common.ncs, backend);        
        } 
    }
    
    START_TIMING;
    // use PTC to solve steady problem
    solv.PseudoTransientContinuation(disc, prec, backend);           
    END_TIMING_DISC(97);    
        
#ifdef TIMING         
    if (disc.common.mpiRank==0) {
        printf("\nComputing initial guess time: %g miliseconds\n", disc.common.timing[99]);   
        printf("Computing AV fields time: %g miliseconds\n", disc.common.timing[98]);   
        printf("Nonlinear solver time: %g miliseconds\n", disc.common.timing[97]);                
    }
#endif    
}

void CSolution::TimeStepping(dstype time, Int istep, Int backend)
{           
    INIT_TIMING;        
           
    // current timestep        
    disc.common.currentstep = istep;
        
    
    for (Int j=0; j<disc.common.tstages; j++) {            
        // current timestage
        disc.common.currentstage = j;
        
        // current time
        if (disc.common.temporalScheme==0)
            disc.common.time = time + disc.common.dt[istep]*disc.common.DIRKcoeff_t[j];
        else
            disc.common.time = time + disc.common.dt[istep]*disc.common.BDFcoeff_t[j];

        if (disc.common.mpiRank==0)
            printf("\nTimestep :  %d,  Timestage :  %d,   Time : %g\n",istep+1,j+1,disc.common.time);            
        
#ifdef TIMING    
        disc.common.timing[100] = 0.0; 
        disc.common.timing[101] = 0.0; 
#endif
        
        // get udg and uhat on the time-dependent boundary 
        if (disc.common.ntdepbc>0) {
            if (j==0)
                GetUdgUhatBou(disc.sol.uhbou0, disc.sol.udgbou0, disc.sol, disc.mesh, disc.common, backend);
            if (j==1)
                GetUdgUhatBou(disc.sol.uhbou1, disc.sol.udgbou1, disc.sol, disc.mesh, disc.common, backend);
            if (j==2)
                GetUdgUhatBou(disc.sol.uhbou2, disc.sol.udgbou2, disc.sol, disc.mesh, disc.common, backend);
        }
        
        START_TIMING;
        // update source term             
        UpdateSource(disc.sol, solv.sys, disc.app, disc.res, disc.common, backend);
        END_TIMING_DISC(100);    
                            
        // solve the problem 
        this->SteadyProblem(backend);                             
                
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
       
    // save solutions into binary files
    SaveSolutions(disc.sol, solv.sys, disc.common, backend);            
}

void CSolution::UnsteadyProblem(Int backend)
{    
    Int nc = disc.common.nc; // number of compoments of (u, q, p)
    Int ncu = disc.common.ncu; // number of compoments of u
    Int ncs = disc.common.ncs;// number of compoments of (s)    
    Int npe = disc.common.npe; // number of nodes on master element    
    Int ne = disc.common.ne1; // number of elements in this subdomain 
    // initial time
    dstype scalar, time = disc.common.time;           
    Int istep;    
        
    if (disc.common.temporalScheme==1)  // BDF schemes
    {        
        // the initial solution at t=0       
        ArrayExtract(solv.sys.udgprev1, disc.sol.udg, npe, nc, ne, 0, npe, 0, ncs, 0, ne, backend);                                                                  

        //DIRK coefficients
        disc.common.temporalScheme = 0; 
        disc.common.tstages = disc.common.torder;                
        TimestepCoefficents(disc.common); 
        
        // Use DIRK scheme to compute the previous solutions to start BDF scheme 
        for (istep=0; istep<disc.common.torder-1; istep++)                  
        {            
            // store previous solutions to calculate the source term        
            PreviousSolutions(disc.sol, solv.sys, disc.common, backend);
            
            // compute the solution at the next step
            this->TimeStepping(time, istep, backend);
                        
            // update time
            time = time + disc.common.dt[istep];        
            
            // extract the current stage solution 
            if (istep==0)
                ArrayExtract(solv.sys.udgprev2, disc.sol.udg, npe, nc, ne, 0, npe, 0, ncs, 0, ne, backend);                                                  
            
            if (istep==1)
                ArrayExtract(solv.sys.udgprev3, disc.sol.udg, npe, nc, ne, 0, npe, 0, ncs, 0, ne, backend);                                                              
        }
         
        //BDF coefficients 
        disc.common.temporalScheme = 1; 
        TimestepCoefficents(disc.common);    
        
        // Now use BDF scheme for the remaining timesteps
        for (istep=disc.common.torder-1; istep<disc.common.tsteps; istep++)            
        {
            // compute the solution at the next step
            this->TimeStepping(time, istep, backend);
            
            // store previous solutions to calculate the source term        
            PreviousSolutions(disc.sol, solv.sys, disc.common, backend);

            // update time
            time = time + disc.common.dt[istep];                    
        }        
    }
    else { // DIRK schemes
        //DIRK coefficients 
        disc.common.temporalScheme = 0; 
        TimestepCoefficents(disc.common); 
                
#ifdef TIMESTEP                  
        struct timeval tv1, tv2;
#endif
        
        Int N = disc.common.ndof1;
        dstype *udgcurt;
        if (disc.common.wave==1)
            udgcurt = &disc.sol.udg[0];
        else
            udgcurt = &solv.sys.u[0];
        ArraySetValue(disc.sol.udgavg, zero, disc.common.ndofsdg1, backend);  
        
        // time stepping with DIRK schemes
        for (istep=0; istep<disc.common.tsteps; istep++)            
        {            
            // store previous solutions to calculate the source term        
            PreviousSolutions(disc.sol, solv.sys, disc.common, backend);

#ifdef TIMESTEP              
            gettimeofday(&tv1, NULL); 
#endif
            
            // compute the solution at the next step
            this->TimeStepping(time, istep, backend);
            
#ifdef TIMESTEP                          
            gettimeofday(&tv2, NULL);            
            if (disc.common.mpiRank==0)
                printf("\nExecution time (in millisec) for timestep %d:  %g\n", istep+1,
                    (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
                    (double)(tv2.tv_sec -tv1.tv_sec )*1000);
#endif
                    
            // compute time-average solution
            ArrayAdd3Vectors(disc.sol.udgavg, solv.sys.udgprev, udgcurt, disc.sol.udgavg, 
                    (dstype) 0.5, (dstype) 0.5, one, disc.common.ndofsdg1, backend);                
                       
//             // update the reduced basis
//             if ((istep==0) && (disc.common.RBdim>0)){                
//                 ArrayExtract(&prec.precond.W[0], udgcurt, npe, ncs, ne, 0, npe, 0, ncu, 0, ne, backend);                 
//                 disc.common.RBcurrentdim = 1;  
//                 PDOT(disc.common.cublasHandle, N, &prec.precond.W[0], inc1, &prec.precond.W[0], inc1, &scalar, backend);     
//                 scalar = sqrt(scalar);
//                 ArrayMultiplyScalar(&prec.precond.W[0], one/scalar, N, backend);                                
//             }
//             else if ((istep==1) && (disc.common.RBdim>1)) {
//                 ArrayExtract(&prec.precond.W[0], udgcurt, npe, ncs, ne, 0, npe, 0, ncu, 0, ne, backend);                                         
//                 ArrayExtract(&prec.precond.W[N], solv.sys.udgprev, npe, ncs, ne, 0, npe, 0, ncu, 0, ne, backend);                 
//                 disc.common.RBcurrentdim = 2;
//                 
//                 PDOT(disc.common.cublasHandle, N, &prec.precond.W[0], inc1, &prec.precond.W[0], inc1, &scalar, backend);     
//                 scalar = sqrt(scalar);
//                 ArrayMultiplyScalar(&prec.precond.W[0], one/scalar, N, backend);                
//                 
//                 PDOT(disc.common.cublasHandle, N, &prec.precond.W[0], inc1, &prec.precond.W[N], inc1, &scalar, backend);     
//                 ArrayAXPBY(&prec.precond.W[N], &prec.precond.W[N], &prec.precond.W[0], one, -scalar, N, backend);                  
//                 PDOT(disc.common.cublasHandle, N, &prec.precond.W[N], inc1, &prec.precond.W[N], inc1, &scalar, backend);        
//                 scalar = sqrt(scalar);
//                 ArrayMultiplyScalar(&prec.precond.W[N], one/scalar, N, backend);                
//             }                
//             else if ((istep>1) && (disc.common.RBdim>2)) {
//                 ArrayExtract(&prec.precond.W[0], udgcurt, npe, ncs, ne, 0, npe, 0, ncu, 0, ne, backend);                                         
//                 ArrayExtract(&prec.precond.W[N], solv.sys.udgprev, npe, ncs, ne, 0, npe, 0, ncu, 0, ne, backend);                 
//                 ArrayExtract(&prec.precond.W[2*N], disc.sol.udgavg, npe, ncs, ne, 0, npe, 0, ncu, 0, ne, backend);         
//                 disc.common.RBcurrentdim = 3;
//                 
//                 PDOT(disc.common.cublasHandle, N, &prec.precond.W[0], inc1, &prec.precond.W[0], inc1, &scalar, backend);     
//                 scalar = sqrt(scalar);
//                 ArrayMultiplyScalar(&prec.precond.W[0], one/scalar, N, backend);                
//                 
//                 PDOT(disc.common.cublasHandle, N, &prec.precond.W[0], inc1, &prec.precond.W[N], inc1, &scalar, backend);     
//                 ArrayAXPBY(&prec.precond.W[N], &prec.precond.W[N], &prec.precond.W[0], one, -scalar, N, backend);                  
//                 PDOT(disc.common.cublasHandle, N, &prec.precond.W[N], inc1, &prec.precond.W[N], inc1, &scalar, backend);        
//                 scalar = sqrt(scalar);
//                 ArrayMultiplyScalar(&prec.precond.W[N], one/scalar, N, backend);  
//                 
//                 PDOT(disc.common.cublasHandle, N, &prec.precond.W[0], inc1, &prec.precond.W[2*N], inc1, &scalar, backend);     
//                 ArrayAXPBY(&prec.precond.W[2*N], &prec.precond.W[2*N], &prec.precond.W[0], one, -scalar, N, backend);                  
//                 PDOT(disc.common.cublasHandle, N, &prec.precond.W[N], inc1, &prec.precond.W[2*N], inc1, &scalar, backend);     
//                 ArrayAXPBY(&prec.precond.W[2*N], &prec.precond.W[2*N], &prec.precond.W[N], one, -scalar, N, backend);                                  
//                 PDOT(disc.common.cublasHandle, N, &prec.precond.W[2*N], inc1, &prec.precond.W[2*N], inc1, &scalar, backend);        
//                 scalar = sqrt(scalar);
//                 ArrayMultiplyScalar(&prec.precond.W[2*N], one/scalar, N, backend);                  
//             }
                
            // update time
            time = time + disc.common.dt[istep];                    
        }           
    }        
}

void CSolution::SolveProblem(Int backend) 
{    
    // compute the geometry quantities
    disc.compGeometry(backend);
    
    // compute the inverse of the mass matrix
    disc.compMassInverse(backend);
        
    // compute q
    if ((disc.common.ncq>0) & (disc.common.wave==0))
        disc.evalQSer(backend);
        
    // set pointer depending on the matrix type
    if (disc.common.precMatrixType==0)
        prec.precond.Cmat = &prec.precond.C[0];
    else //if (disc.common.precMatrixType==2)
        prec.precond.Cmat = &disc.res.Minv[0];    
    
    if (disc.common.PTCparam>0) {
        // set pointer depending on the matrix type
        if (disc.common.ptcMatrixType==0)
            solv.sys.PTCmat = &solv.sys.PTCmatrix[0];
        else //if (disc.common.ptcMatrixType==1)
            solv.sys.PTCmat = &disc.res.Mass[0];        
    }
    
    // copy sol.udg to sys.u
    ArrayExtract(solv.sys.u, disc.sol.udg, disc.common.npe, disc.common.nc, disc.common.ne1, 
            0, disc.common.npe, 0, disc.common.ncu, 0, disc.common.ne1, backend);                                                  
    
    if (disc.common.wave==1)
        // copy sol.wdg to sys.w
        ArrayCopy(solv.sys.w, disc.sol.wdg, disc.common.ndof1, backend);
            
    // initialize uh
    InitUhat(disc.sol, disc.mesh, disc.tmp, disc.common, backend);
        
    if (disc.common.tdep==1) {        
        // allocate memory for uh and udg for the time-dependent boundary conditions
        if (disc.common.ntdepbc>0)
            SetUdgUhatBou(disc.sol, disc.common, backend);                
        
        // solve time-dependent problems
        this->UnsteadyProblem(backend);            
    }
    else {
        // solve steady-state problems
        this->SteadyProblem(backend);        
        
        // save solutions into binary files
        SaveSolutions(disc.sol, solv.sys, disc.common, backend);    
    }    
}

#endif        
