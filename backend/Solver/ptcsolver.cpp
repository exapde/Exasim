/*
    ptcsolver.cpp

    This file contains implementations for nonlinear and linear solvers using Pseudo-Time Continuation (PTC) and Newton's method, 
    with support for reduced basis (RB) preconditioning and timing analysis. The solvers are designed for high-performance 
    scientific computing, potentially leveraging CUDA and MPI for parallelism.

    Functions:

    - int LinearSolver(sysstruct &sys, CDiscretization& disc, CPreconditioner& prec, ofstream &out, Int it, Int backend)
        Solves the linear system arising in each nonlinear iteration. Evaluates the residual, constructs the preconditioner, 
        applies GMRES, and manages timing and logging. Handles reduced basis updates and checks for convergence.

    - void UpdateRB(sysstruct &sys, CDiscretization& disc, CPreconditioner& prec, Int backend)
        Updates the reduced basis vectors used for preconditioning if the current solution increment is significant.

    - Int PTCsolver(sysstruct &sys,  CDiscretization& disc, CPreconditioner& prec, ofstream &out, Int backend)
        Implements the Pseudo-Time Continuation (PTC) nonlinear solver. Iteratively solves the nonlinear system R(u) = 0 
        using the linear solver, updates the solution, checks for divergence, logs residual norms, and manages reduced basis updates.

    - void UpdateRB(sysstruct &sys, CDiscretization& disc, CPreconditioner& prec, Int N, Int backend)
        Alternate version of UpdateRB with explicit dimension argument. Updates the reduced basis vectors for preconditioning.

    - void LinearSolver(sysstruct &sys, CDiscretization& disc, CPreconditioner& prec, ofstream &out, Int N, Int spatialScheme, Int it, Int backend)
        Overloaded LinearSolver supporting different spatial discretization schemes. Assembles the linear system, constructs 
        the preconditioner, applies GMRES, and logs timing information.

    - Int NonlinearSolver(sysstruct &sys,  CDiscretization& disc, CPreconditioner& prec, ofstream &out, Int N, Int spatialScheme, Int backend)
        Implements Newton's method for nonlinear systems. Iteratively solves the linearized system, applies solution updates 
        with damping, computes residuals, manages reduced basis updates, and checks for convergence.

    Key Concepts:
    - Residual evaluation and norm computation for convergence checks.
    - Reduced basis (RB) preconditioning for accelerating linear solves.
    - GMRES iterative solver for linear systems.
    - Timing and performance analysis for solver components.
    - Support for different spatial discretization schemes (e.g., HDG).
    - Logging and debugging support for solution and residual vectors.
    - Handling of divergence and NaN residuals with output and termination.

    Dependencies:
    - CUDA and MPI support (conditional compilation).
    - External array and linear algebra utilities (e.g., ArrayAXPY, ArrayCopy, PNORM).
    - Data structures: sysstruct, CDiscretization, CPreconditioner.

*/
#ifndef __PTCSOLVER
#define __PTCSOLVER

int LinearSolver(sysstruct &sys, CDiscretization& disc, CPreconditioner& prec, ofstream &out, Int it, Int backend)
{    
        
#ifdef TIMING    
    auto begin = chrono::high_resolution_clock::now();      
#endif
   
//#ifdef HAVE_CUDA
//    cudaDeviceSynchronize();
//#endif
     
    // evaluate the residual R(u) and set it to sys.b
    disc.evalResidual(sys.b, sys.u, backend);

    int N = disc.common.ndof1;

    // residual norm
    dstype oldnrm = PNORM(disc.common.cublasHandle, N, sys.b, backend); 
                
    if (disc.common.mpiRank==0 && disc.common.saveResNorm==1 && it==1) {
        disc.common.timing[120] = disc.common.currentstep + 1.0;
        disc.common.timing[121] = disc.common.currentstage + 1.0;
        disc.common.timing[122] = 0.0; 
        disc.common.timing[123] = oldnrm;        
        writearray(out, &disc.common.timing[120], 4);    
    }
    
    // construct the preconditioner
    if (disc.common.RBcurrentdim>0) {
        //prec.ConstructPreconditioner(sys, disc, backend);                  
        prec.ComputeInitialGuessAndPreconditioner(sys, disc, backend); 
        
        // v = u + x 
        //int N = disc.common.ndof1;
        ArrayAXPBY(disc.common.cublasHandle, sys.v, sys.u, sys.x, one, one, N, backend);  
        disc.evalResidual(sys.r, sys.v, backend);  
        dstype nrmr = PNORM(disc.common.cublasHandle, N, sys.r, backend);
                
        if (nrmr>1.05*oldnrm) {
            //ArraySetValue(sys.x, zero, N, backend);
            ArrayMultiplyScalar(disc.common.cublasHandle, sys.x, zero, N, backend);                       
            // reset the reduced basis
            disc.common.RBremovedind = 0;
            disc.common.RBcurrentdim = 0;
            //ArrayCopy(&prec.precond.W[disc.common.RBremovedind*N], sys.x, N, backend);         
            //disc.common.RBcurrentdim = 1;
            //disc.common.RBremovedind = 1;            
        }

        if (nrmr < disc.common.nonlinearSolverTol) 
            return 1;               
    }    
    else {
        //ArraySetValue(sys.x, zero, disc.common.ndof1, backend);
        ArrayMultiplyScalar(disc.common.cublasHandle, sys.x, zero, N, backend);   
    }
    
//   ArraySetValue(sys.x, zero, disc.common.ndof1, backend);
    
#ifdef TIMING        
    auto end = chrono::high_resolution_clock::now();
    disc.common.timing[4] = chrono::duration_cast<chrono::nanoseconds>(end-begin).count()/1e6;        
#endif
    
    // set Wcurrentdim
    //disc.common.Wcurrentdim = disc.common.RBcurrentdim;
    
#ifdef TIMING    
    begin = chrono::high_resolution_clock::now();      
#endif    
    disc.common.linearSolverIter = GMRES(sys, disc, prec, backend);  
    if (disc.common.mpiRank==0)             
        printf("GMRES converges to the tolerance %g within % d iterations and %d RB dimensions\n",disc.common.linearSolverTol,disc.common.linearSolverIter,disc.common.RBcurrentdim);        
    
#ifdef TIMING        
    end = chrono::high_resolution_clock::now();
    disc.common.timing[5] = chrono::duration_cast<chrono::nanoseconds>(end-begin).count()/1e6;        
        
    if (disc.common.mpiRank==0) {
    printf("--------- Linear Solver Analysis -------\n");    
    printf("Constructing preconditioner time: %g miliseconds\n", disc.common.timing[4]);
    printf("Matrix-vector product time: %g miliseconds\n", disc.common.timing[0]);
    printf("Orthgonalization time: %g miliseconds\n", disc.common.timing[1]);
    printf("Applying preconditioner time: %g miliseconds\n", disc.common.timing[2]);
    printf("Solution update time: %g miliseconds\n", disc.common.timing[3]);
    printf("Linear solver time: %g miliseconds\n\n", disc.common.timing[5]);
        
    printf("--------- Residual Calculation Analysis -------\n");
    printf("Copy to buffsend time: %g miliseconds\n", disc.common.timing[13]);
    printf("Non-blocking send/receive time: %g miliseconds\n", disc.common.timing[6]);
    printf("GetUhat time: %g miliseconds\n", disc.common.timing[7]);
    printf("GetQ (interior elements) time: %g miliseconds\n", disc.common.timing[8]);
    printf("RuElem (interior elements) time: %g miliseconds\n", disc.common.timing[9]);
    printf("MPI_WAITALL time: %g miliseconds\n", disc.common.timing[10]);
    printf("RuElem and GetQ (exterior elements) time: %g miliseconds\n", disc.common.timing[11]);
    printf("RuFace time: %g miliseconds\n", disc.common.timing[12]);    
    printf("Copy from buffrecv time: %g miliseconds\n\n", disc.common.timing[14]);    
    
//     printf("--------- Rq Residual Calculation Analysis -------\n");
//     printf("RqElem time: %g miliseconds\n", disc.common.timing[15]);
//     printf("RqFace time: %g miliseconds\n", disc.common.timing[16]);
//     printf("RqApplyMinv time: %g miliseconds\n", disc.common.timing[17]);
//     printf("RqArrayInsert time: %g miliseconds\n", disc.common.timing[18]);
//     
//     printf("RqFace GetArrayAtIndex xdg time: %g miliseconds\n", disc.common.timing[43]);
//     printf("RqFace Node2Gauss xdg time: %g miliseconds\n", disc.common.timing[19]);
//     printf("RqFace facegeom time: %g miliseconds\n", disc.common.timing[20]);
//     printf("RqFace GetElemNodes uhg time: %g miliseconds\n", disc.common.timing[44]);
//     printf("RqFace Node2Gauss uhg time: %g miliseconds\n", disc.common.timing[21]);
//     printf("RqFace ApplyJacNormal time: %g miliseconds\n", disc.common.timing[22]);
//     printf("RqFace Gauss2Node time: %g miliseconds\n", disc.common.timing[23]);
//     printf("RqFace PutFaceNodes time: %g miliseconds\n", disc.common.timing[24]);
//     
//     printf("RqElem GetArrayAtIndex xdg time: %g miliseconds\n", disc.common.timing[41]);
//     printf("RqElem Node2Gauss xdg time: %g miliseconds\n", disc.common.timing[25]);
//     printf("RqElem elemgeom time: %g miliseconds\n", disc.common.timing[26]);
//     printf("RqElem GetArrayAtIndex udg time: %g miliseconds\n", disc.common.timing[42]);
//     printf("RqElem Node2Gauss udg time: %g miliseconds\n", disc.common.timing[27]);
//     printf("RqElem Gauss2Node time: %g miliseconds\n", disc.common.timing[28]);
//     printf("RqElem PutArrayAtIndex time: %g miliseconds\n\n", disc.common.timing[29]);    
//     
//     printf("--------- Ru Residual Calculation Analysis -------\n");
//     printf("RuFace GetArrayAtIndex xdg time: %g miliseconds\n", disc.common.timing[50]);
//     printf("RuFace Node2Gauss xdg time: %g miliseconds\n", disc.common.timing[30]);    
//     printf("RuFace facegeom time: %g miliseconds\n", disc.common.timing[31]);
//     printf("RuFace GetElemNodes uhg time: %g miliseconds\n", disc.common.timing[51]);
//     printf("RuFace Node2Gauss uhg time: %g miliseconds\n", disc.common.timing[52]);
//     printf("RuFace GetArrayAtIndex udg1 time: %g miliseconds\n", disc.common.timing[53]);
//     printf("RuFace Node2Gauss udg1 time: %g miliseconds\n", disc.common.timing[54]);    
//     printf("RuFace GetArrayAtIndex udg2 time: %g miliseconds\n", disc.common.timing[55]);
//     printf("RuFace Node2Gauss udg2 time: %g miliseconds\n", disc.common.timing[56]);    
//     printf("RuFace GetFaceNodes odg1 time: %g miliseconds\n", disc.common.timing[57]);
//     printf("RuFace Node2Gauss odg1 time: %g miliseconds\n", disc.common.timing[58]);    
//     printf("RuFace GetFaceNodes odg2 time: %g miliseconds\n", disc.common.timing[59]);
//     printf("RuFace Node2Gauss odg2 time: %g miliseconds\n", disc.common.timing[60]);                
//     printf("RuFace FhatDriver time: %g miliseconds\n", disc.common.timing[32]);
//     printf("RuFace FbouDriver time: %g miliseconds\n", disc.common.timing[32]);
//     printf("RuFace ApplyJac1 time: %g miliseconds\n", disc.common.timing[61]);       
//     printf("RuFace Gauss2Node time: %g miliseconds\n", disc.common.timing[34]);
//     printf("RuFace PutFaceNodes time: %g miliseconds\n", disc.common.timing[35]);
//     
//     printf("RuElem GetArrayAtIndex xdg time: %g miliseconds\n", disc.common.timing[45]);
//     printf("RuElem Node2Gauss xdg time: %g miliseconds\n", disc.common.timing[36]);
//     printf("RuElem elemgeom time: %g miliseconds\n", disc.common.timing[37]);
//     printf("RuElem GetArrayAtIndex udg time: %g miliseconds\n", disc.common.timing[46]);
//     printf("RuElem Node2Gauss udg time: %g miliseconds\n", disc.common.timing[47]);    
//     printf("RuElem GetElemNodes odg time: %g miliseconds\n", disc.common.timing[48]);
//     printf("RuElem Node2Gauss odg time: %g miliseconds\n", disc.common.timing[62]);    
//     printf("RuElem GetElemNodes sdg time: %g miliseconds\n", disc.common.timing[63]);
//     printf("RuElem Node2Gauss sdg time: %g miliseconds\n", disc.common.timing[64]);
//     printf("RuElem TdfuncDriver time: %g miliseconds\n", disc.common.timing[65]);
//     printf("RuElem SourceDriver time: %g miliseconds\n", disc.common.timing[66]);
//     printf("RuElem ApplyJac1 time: %g miliseconds\n", disc.common.timing[67]);
//     printf("RuElem Gauss2Node source: %g miliseconds\n", disc.common.timing[38]);
//     printf("RuElem FluxDriver time: %g miliseconds\n", disc.common.timing[39]);
//     printf("RuElem Gauss2Node1 time: %g miliseconds\n", disc.common.timing[40]);    
//     printf("RuElem PutArrayAtIndex time: %g miliseconds\n\n", disc.common.timing[49]); 
    }
#endif
    
    return 0;
    // reset Wcurrentdim
    //disc.common.Wcurrentdim = disc.common.RBcurrentdim;
}

void UpdateRB(sysstruct &sys, CDiscretization& disc, CPreconditioner& prec, Int backend)
{
    Int N = disc.common.ndof1;
                    
    dstype nrmr = PNORM(disc.common.cublasHandle, N, sys.x, backend);
    if (nrmr>zero) {
      // update the reduced basis        
      //ArrayCopy(&prec.precond.W[disc.common.RBremovedind*N], sys.x, N, backend);  
      ArrayCopy(disc.common.cublasHandle, &prec.precond.W[disc.common.RBremovedind*N], sys.x, N, backend);  

      // update the current dimension of the RB dimension
      if (disc.common.RBcurrentdim<disc.common.RBdim) 
          disc.common.RBcurrentdim += 1;                    

      // update the position of the RB vector to be replaced  
      disc.common.RBremovedind += 1;
      if (disc.common.RBremovedind==disc.common.RBdim) 
          disc.common.RBremovedind = 0;                
    }
}

Int PTCsolver(sysstruct &sys,  CDiscretization& disc, CPreconditioner& prec, ofstream &out, Int backend)       
{
    Int N = disc.common.ndof1;     
    Int it = 0, maxit = disc.common.nonlinearSolverMaxIter;  
    dstype nrmr, tol;
    tol = disc.common.nonlinearSolverTol; // tolerance for the residual
    
    nrmr = PNORM(disc.common.cublasHandle, N, sys.u, backend);
    if (disc.common.mpiRank==0)
        cout<<"PTC Iteration: "<<it<<",  Solution Norm: "<<nrmr<<endl;                                                    
    
    // compute both the residual vector and sol.udg  
    disc.evalResidual(sys.r, sys.u, backend);
    nrmr = PNORM(disc.common.cublasHandle, N, sys.r, backend);
    if (disc.common.mpiRank==0)
        cout<<"PTC Iteration: "<<it<<",  Residual Norm: "<<nrmr<<endl;                           
    
    // use PTC to solve the system: R(u) = 0
    for (it=0; it<maxit; it++) {                        
        //nrmrold = nrmr;
        
        // solve the linear system: (lambda*B + J(u))x = -R(u)
        int status;
        status = LinearSolver(sys, disc, prec, out, it, backend);
                                
        // update the solution: u = u + x
        ArrayAXPY(disc.common.cublasHandle, sys.u, sys.x, one, N, backend); 

// #ifdef HAVE_CUDA
//        cudaDeviceSynchronize();
// #endif         
       
        // if (status==1)
        //    return it;
       
        // compute both the residual vector and sol.udg  
        disc.evalResidual(sys.r, sys.u, backend);
        nrmr = PNORM(disc.common.cublasHandle, N, sys.r, backend);

        if (nrmr > 1.0e6) {                        
            string filename = disc.common.fileout + "_np" + NumberToString(disc.common.mpiRank) + ".bin";                    
            if (disc.common.saveSolOpt==0)
                writearray2file(filename, sys.u, disc.common.ndof1, backend);
            else
                writearray2file(filename, disc.sol.udg, disc.common.ndofudg1, backend);       

            if (disc.common.mpiRank==0)
                cout<<"Residual is NaN in file ptcsolver.cpp at line 301"<<endl;                           
            
            #ifdef  HAVE_MPI       
                MPI_Finalize();    
            #endif
            
            exit(1);
            //error("\nResidual norm in nonlinear solver is NaN.\n");
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
            UpdateRB(sys, disc, prec, backend);      
        
        // check convergence
        if (nrmr < tol) {            
            return it;   
        }
    }
    
//     if (disc.common.tdep == 1) {
//         disc.common.tdep = 0;
//         disc.evalResidual(sys.r, sys.u, backend);
//         nrmr = PNORM(disc.common.cublasHandle, N, sys.r, backend);
//         cout<<"PTC Iteration: "<<it<<",  Steady-State Residual Norm: "<<nrmr<<endl;   
//         disc.common.tdep = 1;
//     }
    
    return it;
}

void UpdateRB(sysstruct &sys, CDiscretization& disc, CPreconditioner& prec, Int N, Int backend)
{                       
    dstype nrmr = PNORM(disc.common.cublasHandle, N, sys.x, backend);
    if (nrmr>zero) {
      // update the reduced basis        
      ArrayCopy(&prec.precond.W[disc.common.RBremovedind*N], sys.x, N);  
      //ArrayCopy(disc.common.cublasHandle, &prec.precond.W[disc.common.RBremovedind*N], sys.x, N, backend);  

      // update the current dimension of the RB dimension
      if (disc.common.RBcurrentdim<disc.common.RBdim) 
          disc.common.RBcurrentdim += 1;                    

      // update the position of the RB vector to be replaced  
      disc.common.RBremovedind += 1;
      if (disc.common.RBremovedind==disc.common.RBdim) 
          disc.common.RBremovedind = 0;                
    }
}

void LinearSolver(sysstruct &sys, CDiscretization& disc, CPreconditioner& prec, ofstream &out, Int N, Int spatialScheme, Int it, Int backend)
{            
    // evaluate the residual R(u) and set it to sys.b
    if (spatialScheme==0) {
      disc.evalResidual(sys.b, sys.u, backend);
    }
    else if (spatialScheme==1) {            
      auto begin = chrono::high_resolution_clock::now();   
            
      disc.hdgAssembleLinearSystem(sys.b, backend);
              
      auto end = chrono::high_resolution_clock::now();
      double t1 = chrono::duration_cast<chrono::nanoseconds>(end-begin).count()/1e6;      

      if (disc.common.mpiRank==0) printf("hdgAssembleLinearSystem time: %g miliseconds\n", t1);
      
      if (disc.common.debugMode==1) {
        int n = disc.common.npe*disc.common.ncu;
        int m = disc.common.npf*disc.common.nfe*disc.common.ncu;
        int ne = disc.common.ne1;
        writearray2file(disc.common.fileout + NumberToString(it+1) + "newton_AE.bin", disc.res.H, m*m*ne, backend);
        writearray2file(disc.common.fileout + NumberToString(it+1) + "newton_FE.bin", disc.res.Rh, m*ne, backend);
        writearray2file(disc.common.fileout + NumberToString(it+1) + "newton_DUDG.bin", disc.res.Ru, n*ne, backend);
        writearray2file(disc.common.fileout + NumberToString(it+1) + "newton_DUDG_DUH.bin", disc.res.F, n*m*ne, backend);
      }
    }
    
    // construct the preconditioner
    if (disc.common.RBcurrentdim>0) {
        prec.ComputeInitialGuessAndPreconditioner(sys, disc, N, spatialScheme, backend);         
    }    
    else {
        ArraySetValue(sys.x, zero, N);     
    }
            
    auto begin = chrono::high_resolution_clock::now();   
        
    disc.common.linearSolverIter = GMRES(sys, disc, prec, N, spatialScheme, backend);  

    auto end = chrono::high_resolution_clock::now();
    double t1 = chrono::duration_cast<chrono::nanoseconds>(end-begin).count()/1e6;        
    
    if (disc.common.mpiRank==0)  {
        printf("GMRES time: %g miliseconds\n", t1);
        printf("GMRES(%d) converges to the tolerance %g within % d iterations and %d RB dimensions\n",disc.common.gmresRestart,disc.common.linearSolverTol,disc.common.linearSolverIter,disc.common.RBcurrentdim);                    
    }
    
    //if ((disc.common.linearSolverRelError <= disc.common.linearSolverTol*disc.common.linearSolverTolFactor) && (disc.common.mpiRank==0))             
    //    printf("GMRES(%d) converges to the tolerance %g within % d iterations and %d RB dimensions\n",disc.common.gmresRestart,disc.common.linearSolverTol,disc.common.linearSolverIter,disc.common.RBcurrentdim);                    
}

Int NonlinearSolver(sysstruct &sys,  CDiscretization& disc, CPreconditioner& prec, ofstream &out, Int N, Int spatialScheme, Int backend)       
{
    Int it = 0, maxit = disc.common.nonlinearSolverMaxIter;  
    dstype nrmr, nrm0, tol;
    tol = disc.common.nonlinearSolverTol; // tolerance for the residual
                
    //cout<<disc.common.mpiRank<<"  "<<N<<endl;   
    nrmr = PNORM(disc.common.cublasHandle, N, disc.common.ndofuhatinterface, sys.u, backend);
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
      disc.hdgAssembleResidual(sys.b, backend);
            
      nrmr = PNORM(disc.common.cublasHandle, N, disc.common.ndofuhatinterface, sys.b, backend);       
      nrmr += PNORM(disc.common.cublasHandle, disc.common.npe*disc.common.ncu*disc.common.ne1, disc.res.Ru, backend);                 
      if (disc.common.mpiRank==0)
        cout<<"Newton Iteration: "<<0<<",  Residual Norm: "<<nrmr<<endl;          
    }                
    
    // use PTC to solve the system: R(u) = 0
    for (it=0; it<maxit; it++) {              
                      
        // solve the linear system:  J(u) x = -R(u)        
        LinearSolver(sys, disc, prec, out, N, spatialScheme, it, backend);
                        
//         printArray2D(sys.u, disc.common.npf, disc.common.nf, backend);
//         printArray2D(sys.x, disc.common.npf, disc.common.nf, backend);
        
        sys.alpha = 1.0;        
        // update the solution: u = u + alpha*x
        ArrayAXPY(disc.common.cublasHandle, sys.u, sys.x, sys.alpha, N, backend); 
        
        
        if (spatialScheme == 0) {          
          // compute both the residual vector and sol.udg  
          disc.evalResidual(sys.r, sys.u, backend);          
          nrmr = PNORM(disc.common.cublasHandle, N, sys.r, backend);          
        } 
        else if (spatialScheme == 1) {      
          ArrayCopy(disc.sol.uh, sys.u, N);
          hdgGetDUDG(disc.res.Ru, disc.res.F, sys.x, disc.res.Rq, disc.mesh, disc.common, backend);          
          ArrayCopy(sys.v, disc.res.Ru, disc.common.npe*disc.common.ncu*disc.common.ne1);
          UpdateUDG(disc.sol.udg, disc.res.Ru, sys.alpha, disc.common.npe, disc.common.nc, disc.common.ne1, 0, disc.common.npe, 0, disc.common.ncu, 0, disc.common.ne1);                    
                    
          if (disc.common.debugMode==1) {
            writearray2file(disc.common.fileout + NumberToString(it+1) + "newton_x.bin", sys.x, N, backend);
            writearray2file(disc.common.fileout + NumberToString(it+1) + "newton_u.bin", sys.u, N, backend);
            writearray2file(disc.common.fileout + NumberToString(it+1) + "newton_uh.bin", disc.sol.uh, N, backend);
            writearray2file(disc.common.fileout + NumberToString(it+1) + "newton_udg.bin", disc.sol.udg, disc.common.npe*disc.common.nc*disc.common.ne1, backend);
            error("stop for debugging...");
          }          
                    
          if (disc.common.ncq > 0) hdgGetQ(disc.sol.udg, disc.sol.uh, disc.sol, disc.res, disc.mesh, disc.tmp, disc.common, backend);          
          if (disc.common.ncw > 0) GetW(disc.sol.wdg, disc.sol, disc.tmp, disc.app, disc.common, backend);
                              
          nrm0 = nrmr; // original norm          
          // compute the updated residual norm |[Ru; Rh]|
          disc.hdgAssembleResidual(sys.b, backend);          
          nrmr = PNORM(disc.common.cublasHandle, N, disc.common.ndofuhatinterface, sys.b, backend);           
          nrmr += PNORM(disc.common.cublasHandle, disc.common.npe*disc.common.ncu*disc.common.ne1, disc.res.Ru, backend);   
                                    
          // if (nrmr > nrm0 && nrmr > 1.0e6) {                        
          //   string filename = disc.common.fileout + "_np" + NumberToString(disc.common.mpiRank) + ".bin";                    
          //   writearray2file(filename, disc.sol.udg, disc.common.ndofudg1, backend);       
          //   error("Residual norm increases more than 1e6. Save and exit.");                                    
          // }
            
          // damped Newton loop to determine alpha
          while ((nrmr>nrm0 && sys.alpha > 0.1) || IS_NAN(nrmr)) 
          {
            if (disc.common.mpiRank==0)
              printf("Newton Iteration: %d, Alpha: %g, Original Norm: %g,  Updated Norm: %g\n", it+1, sys.alpha, nrm0, nrmr);
            sys.alpha = sys.alpha/2.0;             
            ArrayAXPY(disc.common.cublasHandle, sys.u, sys.x, -sys.alpha, N, backend); 
            ArrayCopy(disc.sol.uh, sys.u, N);
            UpdateUDG(disc.sol.udg, sys.v, -sys.alpha, disc.common.npe, disc.common.nc, disc.common.ne1, 0, disc.common.npe, 0, disc.common.ncu, 0, disc.common.ne1);                    
            if (disc.common.ncq > 0) hdgGetQ(disc.sol.udg, disc.sol.uh, disc.sol, disc.res, disc.mesh, disc.tmp, disc.common, backend);          
            if (disc.common.ncw > 0) GetW(disc.sol.wdg, disc.sol, disc.tmp, disc.app, disc.common, backend);
            disc.hdgAssembleResidual(sys.b, backend);
            nrmr = PNORM(disc.common.cublasHandle, N, disc.common.ndofuhatinterface, sys.b, backend); 
            nrmr += PNORM(disc.common.cublasHandle, disc.common.npe*disc.common.ncu*disc.common.ne1, disc.res.Ru, backend);                       
          }          
        }

        // update the reduced basis space
        ArrayMultiplyScalar(disc.common.cublasHandle, sys.x, sys.alpha, N, backend);   
                        
        if (disc.common.RBdim > 0) UpdateRB(sys, disc, prec, N, backend);         
                
        if (disc.common.mpiRank==0)
          printf("Newton Iteration: %d, Alpha: %g, Original Norm: %g,  Updated Norm: %g\n", it+1, sys.alpha, nrm0, nrmr);
        
        // check convergence
        if (nrmr < tol) return (it+1);           
    }
    
    return it;
}

#endif


