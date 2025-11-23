/*
    solution.cpp

    This file implements the CSolution class methods for solving PDEs using various numerical schemes.
    The main functionalities include initialization, time-stepping, steady-state solving, saving and reading solutions,
    and handling output for both DG and CG methods. The code supports parallel execution (MPI), adaptive timestepping,
    and artificial viscosity (AV) field computation.

    Main Methods:
    -------------
    - SteadyProblem: Solves steady-state problems, computes AV fields, and prepares solution data for output.
    - InitSolution: Initializes solution variables, sets up geometry, mass matrix, and prepares for time-dependent or steady problems.
    - DIRK: Implements time-stepping using Diagonally Implicit Runge-Kutta (DIRK) schemes, including source term updates and solution saving.
    - SteadyProblem_PTC: Solves steady problems using Pseudo-Transient Continuation (PTC) with adaptive timestep control and convergence monitoring.
    - SolveProblem: Entry point for solving either time-dependent or steady-state problems, including solution initialization and output.
    - SaveSolutions / ReadSolutions: Save and load solution data to/from binary files, supporting both time-dependent and steady-state cases.
    - SaveOutputDG / SaveOutputCG: Save DG and CG output data, including post-processing and conversion between DG and CG representations.
    - SaveSolutionsOnBoundary / SaveNodesOnBoundary: Save solution and node data on domain boundaries for post-processing or coupling.

    Features:
    ---------
    - Supports multiple spatial discretization schemes (HDG, DG, etc.).
    - Handles artificial viscosity computation and smoothing, including parallel communication.
    - Adaptive timestep control for PTC and DIRK schemes.
    - MPI parallelization for distributed memory computation.
    - Flexible output and checkpointing for solutions and boundary data.
    - Modular design with external dependencies for geometry, mass matrix, and residual computation.

    Usage:
    ------
    - Instantiate CSolution and call SolveProblem() to solve a PDE problem.
    - Use SaveSolutions(), SaveOutputDG(), SaveOutputCG() for output and post-processing.
    - Configure via disc.common for problem-specific parameters (timestepping, output frequency, AV options, etc.).

    Note:
    -----
    - Requires external files: solution.h, previoussolutions.cpp, updatesolution.cpp, updatesource.cpp, timestepcoeff.cpp, avsolution.cpp.
    - Some methods rely on external functions for array manipulation, MPI communication, and numerical linear algebra.
    - Timing and debugging output are controlled via preprocessor macros (TIMING, TIMESTEP, HAVE_MPI, etc.).
*/
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

template <typename Model>
Int CSolution<Model>::PTCsolver(ofstream &out, Int backend)       
{
    Int N = disc.common.ndof1;     
    Int it = 0, maxit = disc.common.nonlinearSolverMaxIter;  
    dstype nrmr, tol;
    tol = disc.common.nonlinearSolverTol; // tolerance for the residual
    
    nrmr = PNORM(disc.common.cublasHandle, N, solv.sys.u, backend);
    if (disc.common.mpiRank==0)
        cout<<"PTC Iteration: "<<it<<",  Solution Norm: "<<nrmr<<endl;                                                    
    
    // compute both the residual vector and sol.udg  
    disc.evalResidual(solv.sys.r, solv.sys.u, backend);
    nrmr = PNORM(disc.common.cublasHandle, N, solv.sys.r, backend);
    if (disc.common.mpiRank==0)
        cout<<"PTC Iteration: "<<it<<",  Residual Norm: "<<nrmr<<endl;                           
    
    // use PTC to solve the system: R(u) = 0
    for (it=0; it<maxit; it++) {                        

        // solve the linear system: (lambda*B + J(u))x = -R(u)
        int status;
        status = LinearSolver(solv.sys, disc, prec, out, it, backend);
                                
        // update the solution: u = u + x
        ArrayAXPY(disc.common.cublasHandle, solv.sys.u, solv.sys.x, one, N, backend); 

        // compute both the residual vector and sol.udg  
        disc.evalResidual(solv.sys.r, solv.sys.u, backend);
        nrmr = PNORM(disc.common.cublasHandle, N, solv.sys.r, backend);

        if (nrmr > 1.0e6) {   
            string filename = disc.common.fileout + "_np" + NumberToString(disc.common.mpiRank) + ".bin";                    
            writearray2file(filename, disc.sol.udg, disc.common.ndofudg1, backend);       
            if (outsol.is_open()) { outsol.close(); }
            if (outwdg.is_open()) { outwdg.close(); }
            if (outuhat.is_open()) { outuhat.close(); }
            if (outbouxdg.is_open()) { outbouxdg.close(); }
            if (outboundg.is_open()) { outboundg.close(); }
            if (outbouudg.is_open()) { outbouudg.close(); }
            if (outbouwdg.is_open()) { outbouwdg.close(); }
            if (outbouuhat.is_open()) { outbouuhat.close(); }
            if (outqoi.is_open()) { outqoi.close(); }              
            error("Residual norm increases more than 1e6. Save and exit.");                                                
        }
        
        if (disc.common.mpiRank==0 && disc.common.saveResNorm==1) {
            disc.common.timing[122] = it + 0.0; 
            disc.common.timing[123] = nrmr;        
            writearray(out, &disc.common.timing[120], 4);    
        }
        
        if (disc.common.mpiRank==0)
            cout<<"PTC Iteration: "<<it<<",  Residual Norm: "<<nrmr<<endl;                           
                        
        // update the reduced basis
        if ((status==0) && (disc.common.RBdim > 0)) // fix bug here 
            UpdateRB(solv.sys, disc, prec, backend);      
        
        // check convergence
        if (nrmr < tol) {            
            return it;   
        }
    }
        
    return it;
}
template <typename Model>
Int CSolution<Model>::NewtonSolver(ofstream &out, Int N, Int spatialScheme, Int backend)       
{
    Int it = 0, maxit = disc.common.nonlinearSolverMaxIter;  
    dstype nrmr, nrm0, tol;
    tol = disc.common.nonlinearSolverTol; // tolerance for the residual
                
    nrmr = PNORM(disc.common.cublasHandle, N, disc.common.ndofuhatinterface, solv.sys.u, backend);
//     if (disc.common.mpiProcs>1 && disc.common.spatialScheme==1) {
//       dstype nrm = PNORM(disc.common.cublasHandle, disc.common.ncu*disc.common.npf*disc.common.ninterfacefaces, sys.u, backend);
//       nrmr = sqrt(nrmr*nrmr - 0.5*nrm*nrm);
//     }                
    
    if (disc.common.mpiRank==0)
      cout<<"Newton Iteration: "<<it<<",  Solution Norm: "<<nrmr<<endl;                                                        

    if (disc.common.debugMode==1) {
      writearray2file(disc.common.fileout + NumberToString(it) + "newton_uh.bin", disc.sol.uh, N, backend);
      writearray2file(disc.common.fileout + NumberToString(it) + "newton_udg.bin", disc.sol.udg, disc.common.npe*disc.common.nc*disc.common.ne1, backend);
    }

    if (spatialScheme == 1) { 
            
      if (disc.common.ncq > 0) hdgGetQ(disc.sol.udg, disc.sol.uh, disc.sol, disc.res, disc.mesh, disc.tmp, disc.common, backend);          
              
      // compute the residual vector R = [Ru; Rh]
      disc.hdgAssembleResidual(solv.sys.b, backend);
            
      nrmr = PNORM(disc.common.cublasHandle, N, disc.common.ndofuhatinterface, solv.sys.b, backend);       
      nrmr += PNORM(disc.common.cublasHandle, disc.common.npe*disc.common.ncu*disc.common.ne1, disc.res.Ru, backend);                 
      if (disc.common.mpiRank==0)
        cout<<"Newton Iteration: "<<0<<",  Residual Norm: "<<nrmr<<endl;          
    }                
    
    // use PTC to solve the system: R(u) = 0
    for (it=0; it<maxit; it++) {              
                      
        // solve the linear system:  J(u) x = -R(u)        
        LinearSolver(solv.sys, disc, prec, out, N, spatialScheme, it, backend);
        
        solv.sys.alpha = 1.0;        
        // update the solution: u = u + alpha*x
        ArrayAXPY(disc.common.cublasHandle, solv.sys.u, solv.sys.x, solv.sys.alpha, N, backend); 
                
        if (spatialScheme == 0) {          
          // compute both the residual vector and sol.udg  
          disc.evalResidual(solv.sys.r, solv.sys.u, backend);          
          nrmr = PNORM(disc.common.cublasHandle, N, solv.sys.r, backend);          
        } 
        else if (spatialScheme == 1) {      
          ArrayCopy(disc.sol.uh, solv.sys.u, N);
          hdgGetDUDG(disc.res.Ru, disc.res.F, solv.sys.x, disc.res.Rq, disc.mesh, disc.common, backend);          
          ArrayCopy(solv.sys.v, disc.res.Ru, disc.common.npe*disc.common.ncu*disc.common.ne1);
          UpdateUDG(disc.sol.udg, disc.res.Ru, solv.sys.alpha, disc.common.npe, disc.common.nc, disc.common.ne1, 0, disc.common.npe, 0, disc.common.ncu, 0, disc.common.ne1);                    
                    
          if (disc.common.debugMode==1) {
            writearray2file(disc.common.fileout + NumberToString(it+1) + "newton_x.bin", solv.sys.x, N, backend);
            writearray2file(disc.common.fileout + NumberToString(it+1) + "newton_u.bin", solv.sys.u, N, backend);
            writearray2file(disc.common.fileout + NumberToString(it+1) + "newton_uh.bin", disc.sol.uh, N, backend);
            writearray2file(disc.common.fileout + NumberToString(it+1) + "newton_udg.bin", disc.sol.udg, disc.common.npe*disc.common.nc*disc.common.ne1, backend);
            error("stop for debugging...");
          }          
                    
          if (disc.common.ncq > 0) hdgGetQ(disc.sol.udg, disc.sol.uh, disc.sol, disc.res, disc.mesh, disc.tmp, disc.common, backend);          
          if (disc.common.ncw > 0) GetW(disc.sol.wdg, disc.sol, disc.tmp, disc.app, disc.common, backend);
                              
          nrm0 = nrmr; // original norm          
          // compute the updated residual norm |[Ru; Rh]|
          disc.hdgAssembleResidual(solv.sys.b, backend);          
          nrmr = PNORM(disc.common.cublasHandle, N, disc.common.ndofuhatinterface, solv.sys.b, backend);           
          nrmr += PNORM(disc.common.cublasHandle, disc.common.npe*disc.common.ncu*disc.common.ne1, disc.res.Ru, backend);   
                                    
          if (nrmr > nrm0 && nrmr > 1.0e6) {                        
            string filename = disc.common.fileout + "_np" + NumberToString(disc.common.mpiRank) + ".bin";                    
            writearray2file(filename, disc.sol.udg, disc.common.ndofudg1, backend);       
            if (outsol.is_open()) { outsol.close(); }
            if (outwdg.is_open()) { outwdg.close(); }
            if (outuhat.is_open()) { outuhat.close(); }
            if (outbouxdg.is_open()) { outbouxdg.close(); }
            if (outboundg.is_open()) { outboundg.close(); }
            if (outbouudg.is_open()) { outbouudg.close(); }
            if (outbouwdg.is_open()) { outbouwdg.close(); }
            if (outbouuhat.is_open()) { outbouuhat.close(); }
            if (outqoi.is_open()) { outqoi.close(); }              
            error("Residual norm increases more than 1e6. Save and exit.");                                    
          }
            
          // damped Newton loop to determine alpha
          while ((nrmr>nrm0 && solv.sys.alpha > 0.1) || IS_NAN(nrmr)) 
          {
            if (disc.common.mpiRank==0)
              printf("Newton Iteration: %d, Alpha: %g, Original Norm: %g,  Updated Norm: %g\n", it+1, solv.sys.alpha, nrm0, nrmr);
            solv.sys.alpha = solv.sys.alpha/2.0;             
            ArrayAXPY(disc.common.cublasHandle, solv.sys.u, solv.sys.x, -solv.sys.alpha, N, backend); 
            ArrayCopy(disc.sol.uh, solv.sys.u, N);
            UpdateUDG(disc.sol.udg, solv.sys.v, -solv.sys.alpha, disc.common.npe, disc.common.nc, disc.common.ne1, 0, disc.common.npe, 0, disc.common.ncu, 0, disc.common.ne1);                    
            if (disc.common.ncq > 0) hdgGetQ(disc.sol.udg, disc.sol.uh, disc.sol, disc.res, disc.mesh, disc.tmp, disc.common, backend);          
            if (disc.common.ncw > 0) GetW(disc.sol.wdg, disc.sol, disc.tmp, disc.app, disc.common, backend);
            disc.hdgAssembleResidual(solv.sys.b, backend);
            nrmr = PNORM(disc.common.cublasHandle, N, disc.common.ndofuhatinterface, solv.sys.b, backend); 
            nrmr += PNORM(disc.common.cublasHandle, disc.common.npe*disc.common.ncu*disc.common.ne1, disc.res.Ru, backend);                       
          }          
        }

        // update the reduced basis space
        ArrayMultiplyScalar(disc.common.cublasHandle, solv.sys.x, solv.sys.alpha, N, backend);   
                        
        if (disc.common.RBdim > 0) UpdateRB(solv.sys, disc, prec, N, backend);         
                
        if (disc.common.mpiRank==0)
          printf("Newton Iteration: %d, Alpha: %g, Original Norm: %g,  Updated Norm: %g\n", it+1, solv.sys.alpha, nrm0, nrmr);
        
        // check convergence
        if (nrmr < tol) return (it+1);           
    }
    
    return it;
}

template <typename Model>
void CSolution<Model>::SteadyProblem(ofstream &out, Int backend) 
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

        Int nco = disc.common.nco;
        Int ncAV = disc.common.ncAV;
        Int npe = disc.common.npe;
        Int ne = disc.common.ne;            
        
        // store AV field 
        dstype *avField = &disc.res.Rq[0];
        dstype *utm = &disc.res.Rq[npe*ncAV*ne];

        // evaluate AV field
        disc.evalAVfield(avField, backend);

        for (Int iav = 0; iav<disc.common.AVsmoothingIter; iav++){
            // printf("Solution AV smoothing iter: %i\n", iav);
            disc.DG2CG2(avField, avField, utm, disc.common.ncAV, disc.common.ncAV, disc.common.ncAV, backend);

#ifdef  HAVE_MPI    
            Int bsz = disc.common.npe*disc.common.ncAV;
            Int nudg = disc.common.npe*disc.common.nco;
            Int n;

            //for (n=0; n<disc.common.nelemsend; n++)
            //    ArrayCopy(&disc.tmp.buffsend[bsz*n], &disc.sol.odg[nudg*disc.common.elemsend[n]], bsz, backend);
            GetArrayAtIndex(disc.tmp.buffsend, avField, disc.mesh.elemsendodg, bsz*disc.common.nelemsend);

#ifdef HAVE_CUDA
            cudaDeviceSynchronize();
#endif

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
            PutArrayAtIndex(avField, disc.tmp.buffrecv, disc.mesh.elemrecvodg, bsz*disc.common.nelemrecv);
#endif
    //    END_TIMING_DISC(98);    
        }

        // insert avField into odg
        ArrayInsert(disc.sol.odg, avField, npe, nco, ne, 0, npe, nco-ncAV, nco, 0, ne);          
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
      this->PTCsolver(out, backend);           
    }
    else if (disc.common.spatialScheme==1) {
      
      this->NewtonSolver(out, disc.common.ndofuhat, disc.common.spatialScheme, backend);           
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

template <typename Model>
void CSolution<Model>::InitSolution(Int backend) 
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

template <typename Model>
void CSolution<Model>::DIRK(ofstream &out, Int backend)
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
        this->SaveQoI(backend);
        if (vis.savemode > 0) this->SaveParaview(backend); 
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

template <typename Model>
void CSolution<Model>::SteadyProblem_PTC(ofstream &out, Int backend) {

    // initial time
    double time = disc.common.time;           
    double monitor_diff, monitor_scale, delta_monitor;
    int N = disc.common.ndofuhat;
    int NLiters = disc.common.nonlinearSolverMaxIter;
    double nrmr = 0;
    int conv_flag = 0;

    Int nc = disc.common.nc; // number of compoments of (u, q, p)
    Int ncu = disc.common.ncu;// number of compoments of (u)    
    Int ncs = disc.common.ncs;// number of compoments of (s)        
    Int npe = disc.common.npe; // number of nodes on master element    
    //Int ne = common.ne1; // number of elements in this subdomain         
    Int ne2 = disc.common.ne2; // number of elements in this subdomain       
    //Int N = common.ndof1;
    Int N2 = npe*disc.common.ncw*ne2;  

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
                ArrayExtract(disc.sol.udg, solv.sys.udgprev, npe, ncu, ne2, 0, npe, 0, nc, 0, ne2);     

                // Compute a new UH
                GetFaceNodes(disc.sol.uh, disc.sol.udg, disc.mesh.f2e, disc.mesh.perm, disc.common.npf, disc.common.ncu, disc.common.npe, disc.common.nc, disc.common.nf);

                // Recompute gradient from udg old
                hdgGetQ(disc.sol.udg, disc.sol.uh, disc.sol, disc.res, disc.mesh, disc.tmp, disc.common, backend);

                // decrease timestep by 10
                std::cout << "Current time step: " << disc.common.dt[istep] << std::endl;
                disc.common.dt[istep+1] = disc.common.dt[istep] / 10;
                std::cout << "next time step: " << disc.common.dt[istep+1] << std::endl;
                if (disc.common.dt[istep+1] < 1e-8){
                    std::cout << "WARNING: PTC stalled" << std::endl;
                    istep = disc.common.tsteps+1;
                }
            }
            else if (delta_monitor < 0.1 && solv.sys.alpha == 1.0)
            {
                if (disc.common.linearSolverIter < disc.common.linearSolverMaxIter){
                    // increase timestep by 2
                    disc.common.dt[istep+1] = disc.common.dt[istep]*2;
                    // std::cout << "Doubling timestep: " << disc.common.dt[istep+1] << std::endl;
                }
                else{ //TODO: Probably overly conservative, consider turning off 
                    disc.common.dt[istep+1] = disc.common.dt[istep]*1;
                    std::cout << "Too many GMRES iterations, not increasing timestep: " << disc.common.dt[istep+1] << std::endl;
                }
                

                if (delta_monitor < 1e-3 && disc.common.dt[istep] > 1e-4) {
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
                            this->SaveQoI(backend);
                            if (vis.savemode > 0) this->SaveParaview(backend); 
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
            else
            {
                // Do not change  timestep
                disc.common.dt[istep+1] = disc.common.dt[istep];
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
        this->SaveQoI(backend);
        if (vis.savemode > 0) this->SaveParaview(backend); 
        this->SaveSolutionsOnBoundary(backend); 
    }
}

template <typename Model>
void CSolution<Model>::SolveProblem(ofstream &out, Int backend) 
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
        this->SaveQoI(backend);
        if (vis.savemode > 0) this->SaveParaview(backend); 
                
        this->SaveSolutionsOnBoundary(backend);         
        if (disc.common.nce>0)
            this->SaveOutputCG(backend);            
    }        
}

template <typename Model>
void CSolution<Model>::SaveSolutions(Int backend) 
{
    bool save = false;
    if (disc.common.tdep==0) save = true;
    else 
        if (((disc.common.currentstep+1) % disc.common.saveSolFreq) == 0) save = true;             

    if (save == true) {        
        if (disc.common.saveSolOpt==0) {
            if (disc.common.spatialScheme > 0) {
                ArrayExtract(disc.res.Rq, disc.sol.udg, disc.common.npe, disc.common.nc, disc.common.ne1, 0, disc.common.npe, 0, disc.common.ncu, 0, disc.common.ne1);                                                  
                writearray(outsol, disc.res.Rq, disc.common.ndof1, backend);    
            }
            else
                writearray(outsol, solv.sys.u, disc.common.ndof1, backend);    
        }
        else
            writearray(outsol, disc.sol.udg, disc.common.ndofudg1, backend);    
        
        if (disc.common.ncw>0)
            writearray(outwdg, disc.sol.wdg, disc.common.ndofw1, backend);

        if (disc.common.spatialScheme==1)
            writearray(outuhat, disc.sol.uh, disc.common.ndofuhat, backend);

        if (disc.common.compudgavg == 1) {
            string fn1 = disc.common.fileout + "solavg_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin"; 
            writearray2file(fn1, disc.sol.udgavg, disc.common.ndofudg1+1, backend);
        }        
    }

    

   // if (disc.common.tdep==1) { 
   //      if (((disc.common.currentstep+1) % disc.common.saveSolFreq) == 0)             
   //      {        
   //          string filename = disc.common.fileout + "udg_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";     
   //          if (disc.common.saveSolOpt==0)
   //              writearray2file(filename, solv.sys.u, disc.common.ndof1, backend);
   //          else
   //              writearray2file(filename, disc.sol.udg, disc.common.ndofudg1, backend);
   // 
   //          if (disc.common.ncw>0) {
   //              string fn = disc.common.fileout + "_wdg_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";                    
   //              writearray2file(fn, solv.sys.wtmp, disc.common.ndofw1, backend);
   //          }                        
   // 
   //          if (disc.common.compudgavg == 1) {
   //              string fn1 = disc.common.fileout + "avg_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin"; 
   //              writearray2file(fn1, disc.sol.udgavg, disc.common.ndofudg1+1, backend);
   //          }
   // 
   //          if (disc.common.spatialScheme==1) {
   //              string fn2 = disc.common.fileout + "_uhat_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";                    
   //              writearray2file(fn2, disc.sol.uh, disc.common.ndofuhat, backend);        
   //          }
   //      }    
   // }
   // else {
   //      string filename = disc.common.fileout + "udg_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";                    
   //      if (disc.common.saveSolOpt==0)
   //          writearray2file(filename, solv.sys.u, disc.common.ndof1, backend);
   //      else
   //          writearray2file(filename, disc.sol.udg, disc.common.ndofudg1, backend);       
   // 
   //      if (disc.common.ncw>0) {
   //          string fn = disc.common.fileout + "_wdg_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";                    
   //          writearray2file(fn, disc.sol.wdg, disc.common.ndofw1, backend);     
   //      }                
   // 
   //      if (disc.common.spatialScheme==1) {
   //          string filename1 = disc.common.fileout + "_uhat_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";                    
   //          writearray2file(filename1, disc.sol.uh, disc.common.ndofuhat, backend);        
   //      }
   // }    
}

template <typename Model>
void CSolution<Model>::ReadSolutions(Int backend) 
{
   if (disc.common.tdep==1) { 
        if (((disc.common.currentstep+1) % disc.common.saveSolFreq) == 0)             
        {        
            string filename = disc.common.fileout + "udg_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";     
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
        string filename = disc.common.fileout + "udg_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";                    
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
 

template <typename Model>
void CSolution<Model>::SaveParaview(Int backend) 
{
    // Decide whether we should write a file on this step
    bool writeSolution = false;
    
    if (disc.common.tdep == 1) {
       if (disc.common.currentstep==0 && disc.common.mpiRank==0) {
          string ext = (disc.common.mpiProcs==1) ? "vtu" : "pvtu";                                  
          vis.pvdwrite_series(disc.common.fileout + "vis", disc.common.dt, disc.common.tsteps, disc.common.saveSolFreq, ext);                          
       }
        
        // Time-dependent: only write every 'saveSolFreq' steps
        writeSolution = ((disc.common.currentstep + 1) % disc.common.saveSolFreq) == 0;
    } else {
        // Steady / not time-dependent: always write
        writeSolution = true;
    }

   if (writeSolution) { 
       int nc = disc.common.nc;  
       int ncx = disc.common.ncx;   
       int nco = disc.common.nco;  
       int ncw = disc.common.ncw;  
       int nsca = disc.common.nsca; 
       int nvec = disc.common.nvec;  
       int nten = disc.common.nten;     
       int npe  = disc.common.npe;     
       int ne   = disc.common.ne1;      
       int ndg  = npe * ne;
       int ncg  = vis.npoints;
    
       dstype* xdg = &disc.tmp.tempn[0];  
       dstype* udg = disc.res.Rq;   
       dstype* vdg = &disc.tmp.tempn[npe*ncx*ne];    
       dstype* wdg = disc.res.Ru;     
       dstype* f = solv.sys.v;
    
       GetElemNodes(xdg, disc.sol.xdg, npe, ncx, 0, ncx, 0, ne);
       GetElemNodes(udg, disc.sol.udg, npe, nc, 0, nc, 0, ne);
       if (nco > 0) GetElemNodes(vdg, disc.sol.odg, npe, nco, 0, nco, 0, ne);
       if (ncw > 0) GetElemNodes(wdg, disc.sol.wdg, npe, ncw, 0, ncw, 0, ne);
    
       if (nsca > 0) {        
            VisScalarsDriver(f, xdg, udg, vdg, wdg, disc.mesh, disc.master, disc.app, disc.sol, disc.tmp, disc.common, npe, 0, ne, backend);                                 
            VisDG2CG(vis.scafields, f, disc.mesh.cgent2dgent, disc.mesh.colent2elem, disc.mesh.rowent2elem, ne, ncg, ndg, 1, 1, nsca);
       }    
       if (nvec > 0) {        
            VisVectorsDriver(f, xdg, udg, vdg, wdg, disc.mesh, disc.master, disc.app, disc.sol, disc.tmp, disc.common, npe, 0, ne, backend);                                 
            VisDG2CG(vis.vecfields, f, disc.mesh.cgent2dgent, disc.mesh.colent2elem, disc.mesh.rowent2elem, ne, ncg, ndg, 3, ncx, nvec);
       }
       if (nten > 0) {        
            VisTensorsDriver(f, xdg, udg, vdg, wdg, disc.mesh, disc.master, disc.app, disc.sol, disc.tmp, disc.common, npe, 0, ne, backend);                                 
            VisDG2CG(vis.tenfields, f, disc.mesh.cgent2dgent, disc.mesh.colent2elem, disc.mesh.rowent2elem, ne, ncg, ndg, vis.ntc, vis.ntc, nvec);
       }

       string baseName = disc.common.fileout + "vis";
       if (disc.common.tdep == 1) {
           std::ostringstream ss; 
           ss << std::setw(6) << std::setfill('0') << disc.common.currentstep+disc.common.timestepOffset+1; 
           baseName = baseName + "_" + ss.str();           
       }
       
       if (disc.common.mpiProcs==1)                
            vis.vtuwrite(baseName, vis.scafields, vis.vecfields, vis.tenfields);
       else
            vis.vtuwrite_parallel(baseName, disc.common.mpiRank, disc.common.mpiProcs, vis.scafields, vis.vecfields, vis.tenfields);
   }
}

template <typename Model>
void CSolution<Model>::SaveQoI(Int backend) 
{
    if (disc.common.nvqoi > 0) qoiElement(disc.sol, disc.res, disc.app, disc.master, disc.mesh, disc.tmp, disc.common);
    if (disc.common.nsurf > 0) qoiFace(disc.sol, disc.res, disc.app, disc.master, disc.mesh, disc.tmp, disc.common);

    if (disc.common.mpiRank==0 && (disc.common.nvqoi > 0 || disc.common.nsurf > 0)) {
        if (disc.common.tdep==1) 
            outqoi << std::setw(16) << std::scientific << std::setprecision(6) << disc.common.time;
        else outqoi << std::setw(16) << std::scientific << std::setprecision(6) << 0.0;
        for (size_t i = 0; i < disc.common.nvqoi; ++i) 
            outqoi << std::setw(16) << std::scientific << std::setprecision(6) << disc.common.qoivolume[i];            
        for (size_t i = 0; i < disc.common.nsurf; ++i) 
            outqoi << std::setw(16) << std::scientific << std::setprecision(6) << disc.common.qoisurface[i];            
        outqoi << "\n";
    }
}

template <typename Model>
void CSolution<Model>::SaveOutputDG(Int backend) 
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

template <typename Model>
void CSolution<Model>::SaveOutputCG(Int backend) 
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

template <typename Model>
void CSolution<Model>::SaveSolutionsOnBoundary(Int backend) 
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
                    Int npe = disc.common.npe; // number of nodes on master face      
                    Int nf = f2-f1;
                    Int nn = npf*nf; 
                    Int nc = disc.common.nc; // number of compoments of (u, q, p)            
                    Int ncu = disc.common.ncu;
                    Int ncw = disc.common.ncw;
                    GetArrayAtIndex(disc.tmp.tempn, disc.sol.udg, &disc.mesh.findudg1[npf*nc*f1], nn*nc);
                    writearray(outbouudg, disc.tmp.tempn, nn*nc, backend);
                    GetElemNodes(disc.tmp.tempn, disc.sol.uh, npf, ncu, 0, ncu, f1, f2);
                    writearray(outbouuhat, disc.tmp.tempn, nn*ncu, backend);
                    if (ncw>0) {
                        GetFaceNodes(disc.tmp.tempn, disc.sol.wdg, disc.mesh.facecon, npf, ncw, npe, ncw, f1, f2, 1);      
                        writearray(outbouwdg, disc.tmp.tempn, nn*ncw, backend);
                    }

                }
            }          
        }                                
    }
}

template <typename Model>
void CSolution<Model>::SaveNodesOnBoundary(Int backend) 
{   
    if ( disc.common.saveSolBouFreq>0 ) {
        for (Int j=0; j<disc.common.nbf; j++) {
            Int f1 = disc.common.fblks[3*j]-1;
            Int f2 = disc.common.fblks[3*j+1];    
            Int ib = disc.common.fblks[3*j+2];            
            if (ib == disc.common.ibs) {     
                Int nd = disc.common.nd; 
                Int npf = disc.common.npf; // number of nodes on master face      
                Int nf = f2-f1;
                Int nn = npf*nf; 
                Int ncx = disc.common.ncx; // number of compoments of (u, q, p)                            
                GetArrayAtIndex(disc.tmp.tempn, disc.sol.xdg, &disc.mesh.findxdg1[npf*ncx*f1], nn*ncx);                
                writearray(outbouxdg, disc.tmp.tempn, nn*ncx, backend);

                Int n1 = nn*ncx;                           // nlg
                Int n2 = nn*(ncx+nd);                      // jac
                Int n3 = nn*(ncx+nd+1);                    // Jg
                if (nd==1) {
                    FaceGeom1D(&disc.tmp.tempn[n2], &disc.tmp.tempn[n1], &disc.tmp.tempn[n3], nn);    
                    FixNormal1D(&disc.tmp.tempn[n1], &disc.mesh.facecon[2*f1], nn);    
                }
                else if (nd==2){
                    Node2Gauss(disc.common.cublasHandle, &disc.tmp.tempn[n3], disc.tmp.tempn, &disc.master.shapfnt[npf*npf], npf, npf, nf*nd, backend);                
                    FaceGeom2D(&disc.tmp.tempn[n2], &disc.tmp.tempn[n1], &disc.tmp.tempn[n3], nn);
                }
                else if (nd==3) {
                    Node2Gauss(disc.common.cublasHandle, &disc.tmp.tempn[n3], disc.tmp.tempn, &disc.master.shapfnt[npf*npf], npf, npf, nf*nd, backend);                     
                    Node2Gauss(disc.common.cublasHandle, &disc.tmp.tempn[n3+nn*nd], disc.tmp.tempn, &disc.master.shapfnt[2*npf*npf], npf, npf, nf*nd, backend);                
                    FaceGeom3D(&disc.tmp.tempn[n2], &disc.tmp.tempn[n1], &disc.tmp.tempn[n3], nn);
                }
                writearray(outboundg, &disc.tmp.tempn[n1], nn*nd, backend);
            }
        }              
    }
}

// void CSolution::SaveSolutionsOnBoundary(Int backend) 
// {   
//     if ( disc.common.saveSolBouFreq>0 ) {
//         if (((disc.common.currentstep+1) % disc.common.saveSolBouFreq) == 0)             
//         {        
//             for (Int j=0; j<disc.common.nbf; j++) {
//                 Int f1 = disc.common.fblks[3*j]-1;
//                 Int f2 = disc.common.fblks[3*j+1];    
//                 Int ib = disc.common.fblks[3*j+2];            
//                 if (ib == disc.common.ibs) {     
//                     Int npf = disc.common.npf; // number of nodes on master face      
//                     Int nf = f2-f1;
//                     Int nn = npf*nf; 
//                     Int nc = disc.common.nc; // number of compoments of (u, q, p)            
//                     GetArrayAtIndex(disc.tmp.tempn, disc.sol.udg, &disc.mesh.findudg1[npf*nc*f1], nn*nc);
//                     string filename = disc.common.fileout + "bou_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";     
//                     writearray2file(filename, disc.tmp.tempn, nn*nc, backend);            
//                 }
//             }                                               
//         }                                
//     }
// }
// 
// void CSolution::SaveNodesOnBoundary(Int backend) 
// {   
//     if ( disc.common.saveSolBouFreq>0 ) {
//         for (Int j=0; j<disc.common.nbf; j++) {
//             Int f1 = disc.common.fblks[3*j]-1;
//             Int f2 = disc.common.fblks[3*j+1];    
//             Int ib = disc.common.fblks[3*j+2];            
//             if (ib == disc.common.ibs) {     
//                 Int npf = disc.common.npf; // number of nodes on master face      
//                 Int nf = f2-f1;
//                 Int nn = npf*nf; 
//                 Int ncx = disc.common.ncx; // number of compoments of (u, q, p)                            
//                 GetArrayAtIndex(disc.tmp.tempn, disc.sol.xdg, &disc.mesh.findxdg1[npf*ncx*f1], nn*ncx);
//                 string filename = disc.common.fileout + "node_np" + NumberToString(disc.common.mpiRank-disc.common.fileoffset) + ".bin";     
//                 writearray2file(filename, disc.tmp.tempn, nn*ncx, backend);            
//             }
//         }                                                                   
//     }
// }

#endif        
