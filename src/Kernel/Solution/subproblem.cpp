//*************** Solver for the PDE model without subproblem **************//

void ApplyPreconditioner(dstype* x, CDiscretization& disc, CPreconditioner& pre, Int backend)
{
    Int N = disc.common.ndof1;            
    // compute  x = C*x
    ArrayCopy(disc.res.Ru, x, N, backend);
    ApplyMatrix(disc.common.cublasHandle, x, pre.precond.Cmat, disc.res.Ru, 
        disc.common.npe, disc.common.ncu, disc.common.ne1, 
        disc.common.precMatrixType, disc.common.curvedMesh, backend);           
}

void Residual(dstype* Ru, dstype* u, CDiscretization& disc, Int backend)
{
    // insert u into udg
    ArrayInsert(disc.sol.udg, u, disc.common.npe, disc.common.nc, disc.common.ne, 0, disc.common.npe, 
            0, disc.common.ncu, 0, disc.common.ne1, backend);  

    // compute the residual vector R(u)
    Residual(disc.sol, disc.res, disc.app, disc.master, disc.mesh, disc.tmp, disc.common, 
            disc.common.cublasHandle, backend);

    // copy the residual vector to Ru
    ArrayCopy(Ru, disc.res.Ru, disc.common.ndof1, backend);          
}

void MatVec(dstype *w, dstype *v, dstype *u, dstype *Ru, CDiscretization& disc, Int backend)
{
  MatVec(w, disc.sol, disc.res, disc.app, disc.master, disc.mesh, disc.tmp, disc.common, 
            disc.common.cublasHandle, v, u, Ru, backend); 
}

void ComputeInitialGuessAndPreconditioner(sysstruct& sys, CDiscretization& disc, 
        CPreconditioner& pre, Int backend)
{           
    Int N = disc.common.ndof1;
    Int RBdim = disc.common.RBcurrentdim;
    dstype *RBcoef = &disc.tmp.tempn[0];
    dstype *RBcoef_tmp = &disc.tmp.tempn[RBdim];
    dstype *H_tmp = &disc.tmp.tempg[RBdim*RBdim];
    
    // compute V = J(u)*W        
    int i = disc.common.RBremovedind -1;
    if (disc.common.RBremovedind==0)
         i = RBdim-1;
    //disc.evalMatVec(&pre.precond.U[i*N], &pre.precond.W[i*N], sys.u, sys.b, backend);      
    MatVec(&pre.precond.U[i*N], &pre.precond.W[i*N], sys.u, sys.b, disc, backend);
    
    /* RBcoef_tmp = V^T b */
    PGEMTV(disc.common.cublasHandle, N, RBdim, &one, pre.precond.U, N, sys.b, inc1, 
            &zero, RBcoef_tmp, inc1, H_tmp, backend);
    
    /* H = V^T V */
    PGEMTM(disc.common.cublasHandle, RBdim, RBdim, N, &one, pre.precond.U, N, pre.precond.U, N, &zero, 
            pre.precond.H, RBdim, H_tmp, backend);    
    
    /* Solve the linear system */
    Inverse(disc.common.cublasHandle, pre.precond.H, H_tmp, pre.precond.ipiv, RBdim, 1, backend);       

    // RBcoef = inverse(H)*RBcoef_tmp
    PGEMNV(disc.common.cublasHandle, RBdim, RBdim, &one, pre.precond.H, RBdim, 
            RBcoef_tmp, inc1, &zero, RBcoef, inc1, backend);    
        
    // compute x = -W*RBcoef
    PGEMNV(disc.common.cublasHandle, N, RBdim, &minusone, pre.precond.W, N, RBcoef, 
            inc1, &zero, sys.x, inc1, backend);                                                 
}

int GMRESsub(sysstruct &sys, CDiscretization &disc, CPreconditioner& prec, Int backend)
{
    Int maxit, nrest, orthogMethod, n1, i, k, j = 0;
    Int ncu = disc.common.ncu;
    Int npe = disc.common.npe;
    Int ne = disc.common.ne1;
    Int N = npe*ncu*ne;
    dstype nrmb, nrmr, tol, scalar;
    tol = min(0.01,disc.common.linearSolverTol*disc.common.linearSolverTolFactor);
    maxit = disc.common.linearSolverMaxIter;
    nrest = disc.common.gmresRestart;
    orthogMethod = disc.common.gmresOrthogMethod;
    n1 = nrest + 1;
                
    dstype *s, *y, *cs, *sn, *H;
    s = &sys.tempmem[0];
    y = &sys.tempmem[n1];
    cs = &sys.tempmem[2*n1];
    sn = &sys.tempmem[3*n1];
    H = &sys.tempmem[4*n1];
    
    // calculate Ritz values of polynomial preconditioner
    if  (disc.common.ppdegree > 1) {
        getPoly(disc, prec, sys, sys.lam, sys.randvect, sys.ipiv, N, disc.common.ppdegree, backend);        
        dstype lmin=1.0e10, lmax=-1.0e10, smax=0.0; 
        int m = disc.common.ppdegree;
        for (int i=0; i<m; i++) {
            lmax = (sys.lam[2*m+i] > lmax) ? sys.lam[2*m+i] : lmax;
            lmin = (sys.lam[2*m+i] < lmin) ? sys.lam[2*m+i] : lmin;
            smax = (fabs(sys.lam[3*m+i]) > smax) ? fabs(sys.lam[3*m+i]) : smax;
        }        
        if (smax == 0.0) {
            for (int i=1; i<=m; i++) 
                sys.lam[i-1] = lmin + (lmax-lmin)*(0.5 - 0.5*cos((2*i-1)*M_PI/(2*m))/cos(M_PI/(2*m)));            
            LejaSort(&sys.lam[2*m], &sys.lam[3*m], sys.lam, &sys.lam[3*m], &sys.lam[4*m], m);
        }        
//         for (int i=0; i<disc.common.ppdegree; i++)
//             if (disc.common.mpiRank==0) cout<<sys.lam[2*disc.common.ppdegree+i]<<"  "<<sys.lam[3*disc.common.ppdegree+i]<<endl;
    }    
    
    // compute b
    nrmb = PNORM(disc.common.cublasHandle, N, sys.b, backend);    
    
    scalar = PNORM(disc.common.cublasHandle, N, sys.x, backend);    
    if (scalar>1e-12) {
        // r = A*x
        //disc.evalMatVec(sys.r, sys.x, sys.u, sys.b, backend);
        MatVec(sys.r, sys.x, sys.u, sys.b, disc, backend);

        // compute the new RHS vector: r = -b - A*x
        ArrayAXPBY(sys.r, sys.b, sys.r, minusone, minusone, N, backend);    

        // norm of the new RHS vector
        nrmr = PNORM(disc.common.cublasHandle, N, sys.r, backend);
    }
    else {
        ArrayAXPBY(sys.r, sys.b, sys.b, minusone, zero, N, backend);
        nrmr = nrmb;  
    }

//     if (disc.common.mpiRank==0)
//         cout<<"Old RHS Norm: "<<nrmb<<",  New RHS Norm: "<<nrmr<<endl; 
    
    disc.common.linearSolverTolFactor = nrmb/nrmr;
    tol = min(0.1,disc.common.linearSolverTol*disc.common.linearSolverTolFactor);
        
    //disc.common.ppdegree = 0;
    
    // compute r = P*r
    if (disc.common.ppdegree>1) {
        ApplyPreconditioner(sys.r, disc, prec, backend);
        //ApplyComponentNorm(disc, sys.normcu, sys.r, disc.res.Ru, ncu, ncu, npe, ne, backend);
        ApplyPoly(sys.r, disc, prec, sys, sys.q, sys.p, N, backend);    
    }
    else if (disc.common.RBcurrentdim>=0) {
        ApplyPreconditioner(sys.r, disc, prec, backend);
        //ApplyComponentNorm(disc, sys.normcu, sys.r, disc.res.Ru, ncu, ncu, npe, ne, backend);
    }
    
    // compute ||r||
    nrmb = PNORM(disc.common.cublasHandle, N, sys.r, backend);
    nrmr = nrmb;
        
    j = 1;
    while (j < maxit) {
        // v = r/||r||
        ArrayAXPB(sys.v, sys.r, one/nrmr, zero, N, backend);
        
        // initialize s
        s[0] = nrmr;
        for (k=1; k<n1; k++)
            s[k] = zero;
        
        for (i = 0; i < nrest && j < maxit; i++, j++) {
            Int m = i + 1;
            
            // compute v[m] = A*v[i]
            START_TIMING;                        
            //disc.evalMatVec(&sys.v[m*N], &sys.v[i*N], sys.u, sys.b, backend);         
            MatVec(&sys.v[m*N], &sys.v[i*N], sys.u, sys.b, disc, backend);
            END_TIMING_DISC(0);    
                                    
            START_TIMING;                                    
            // compute v[m] = P*v[m]
            if (disc.common.ppdegree>1) {
                ApplyPreconditioner(&sys.v[m*N], disc, prec, backend);   
                //ApplyComponentNorm(disc, sys.normcu, &sys.v[m*N], disc.res.Ru, ncu, ncu, npe, ne, backend);
                ApplyPoly(&sys.v[m*N], disc, prec, sys, sys.q, sys.p, N, backend);
            }
            else if (disc.common.RBcurrentdim>=0) {
                ApplyPreconditioner(&sys.v[m*N], disc, prec, backend);   
                //ApplyComponentNorm(disc, sys.normcu, &sys.v[m*N], disc.res.Ru, ncu, ncu, npe, ne, backend);
            }
            END_TIMING_DISC(2);    
                        
            // orthogonalize Krylov vectors
            START_TIMING;    
            if (orthogMethod == 0)
                MGS(disc.common.cublasHandle, sys.v, &H[n1*i], N, m, backend);
            else
                CGS(disc.common.cublasHandle, sys.v, &H[n1*i], y, N, m, backend);
            END_TIMING_DISC(1);    
            
            // Apply Givens rotation to compute s
            ApplyGivensRotation(&H[n1*i], s, cs, sn, i, 0);
            
            // compute relative error
            disc.common.linearSolverRelError = fabs(s[i+1])/nrmb;
            
            // check convergence and update solution: x = x + v*s
            if (disc.common.linearSolverRelError < tol) {                
                UpdateSolution(disc.common.cublasHandle, sys.x, y, H, s, sys.v, i, N, n1, backend);
                return j;
            }            
        }        
        
        // update solution: x = x + v*s
        UpdateSolution(disc.common.cublasHandle, sys.x, y, H, s, sys.v, nrest-1, N, n1, backend);
               
        // compute r = A*x
        //disc.evalMatVec(sys.r, sys.x, sys.u, sys.b, backend);
        MatVec(sys.r, sys.x, sys.u, sys.b, disc, backend);
        
        // r = -b - A*x
        ArrayAXPBY(sys.r, sys.b, sys.r, minusone, minusone, N, backend);
        
        // r = P*r = P*(-b-A*x)
        if (disc.common.ppdegree>1) {
            ApplyPreconditioner(sys.r, disc, prec, backend);            
            //ApplyComponentNorm(disc, sys.normcu, sys.r, disc.res.Ru, ncu, ncu, npe, ne, backend);
            ApplyPoly(sys.r, disc, prec, sys, sys.q, sys.p, N, backend);
        }
        else if (disc.common.RBcurrentdim>=0) {
            ApplyPreconditioner(sys.r, disc, prec, backend);            
            //ApplyComponentNorm(disc, sys.normcu, sys.r, disc.res.Ru, ncu, ncu, npe, ne, backend);
        }
        
        // compute relative error
        nrmr = PNORM(disc.common.cublasHandle, N, sys.r, backend);
        disc.common.linearSolverRelError = nrmr/nrmb;
        
        // check convergence
        if (disc.common.linearSolverRelError < tol) {
            return j;
        }
    }       
    
//     if (disc.common.linearSolverRelError > tol) {
//         if (disc.common.mpiRank==0) {
//             printf("Warning: GMRES(%d) does not converge to the tolerance %g within % d iterations\n",nrest,tol,maxit);
//             printf("Warning: The current relative error is %g \n",disc.common.linearSolverRelError);
//         }
//     }
    return j;
}

int LinearSolver(sysstruct &sys, CDiscretization& disc, CPreconditioner& prec, Int backend)
{            
    // evaluate the residual R(u) and set it to sys.b
    //disc.evalResidual(sys.b, sys.u, backend);
    Residual(sys.b, sys.u, disc, backend);

    int N = disc.common.ndof1;

    // residual norm
    dstype oldnrm = PNORM(disc.common.cublasHandle, N, sys.b, backend); 
                   
    if (oldnrm < disc.common.nonlinearSolverTol) {
      ArraySetValue(sys.x, zero, disc.common.ndof1, backend);
      return 1;           
    }
    
    // construct the preconditioner
    if (disc.common.RBcurrentdim>0) {
        ComputeInitialGuessAndPreconditioner(sys, disc, prec, backend); 
        
        // v = u + x 
        ArrayAXPBY(sys.v, sys.u, sys.x, one, one, N, backend);  
        //disc.evalResidual(sys.r, sys.v, backend);  
        Residual(sys.r, sys.v, disc, backend);
        dstype nrmr = PNORM(disc.common.cublasHandle, N, sys.r, backend);
                
        if (nrmr>1.05*oldnrm) {
            ArraySetValue(sys.x, zero, disc.common.ndof1, backend);
            // reset the reduced basis
            disc.common.RBremovedind = 0;
            disc.common.RBcurrentdim = 0;
        }

        if (nrmr < disc.common.nonlinearSolverTol) 
            return 1;               
    }    
    else
        ArraySetValue(sys.x, zero, disc.common.ndof1, backend);
    
    disc.common.linearSolverIter = GMRESsub(sys, disc, prec, backend);  
    
//     if ((disc.common.linearSolverRelError <= disc.common.linearSolverTol*disc.common.linearSolverTolFactor) && (disc.common.mpiRank==0))             
//         printf("GMRES converges to the tolerance %g within % d iterations and %d RB dimensions\n",disc.common.linearSolverTol,disc.common.linearSolverIter,disc.common.RBcurrentdim);        
                 
    return 0;
}

Int NewtonSolver(sysstruct &sys, CDiscretization& disc, CPreconditioner& prec, Int backend)       
{
    Int N = disc.common.ndof1;     
    Int it = 0, maxit = disc.common.nonlinearSolverMaxIter;  
    dstype nrmr, tol;
    tol = disc.common.nonlinearSolverTol; // tolerance for the residual
        
    nrmr = PNORM(disc.common.cublasHandle, N, sys.u, backend);
//     if (disc.common.mpiRank==0)
//         cout<<"Newton Iteration: "<<it<<",  Solution Norm: "<<nrmr<<endl;                                                    
    
    // use PTC to solve the system: R(u) = 0
    for (it=1; it<maxit; it++) {                        
        //nrmrold = nrmr;
        
        // solve the linear system: (lambda*B + J(u))x = -R(u)
        int status;
        status = LinearSolver(sys, disc, prec, backend);
                                
        // update the solution: u = u + x
        ArrayAXPBY(sys.u, sys.u, sys.x, one, one, N, backend); 

        if (status==1)
           return it;
       
        // compute both the residual vector and sol.udg  
        //disc.evalResidual(sys.r, sys.u, backend);
        Residual(sys.r, sys.u, disc, backend);
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
        }
                
//         if (disc.common.mpiRank==0)
//             cout<<"Newton Iteration: "<<it<<",  Residual Norm: "<<nrmr<<endl;                           
                        
        // update the reduced basis
        if (status==0)
            UpdateRB(sys, disc, prec, backend);      
        
        // check convergence
        if (nrmr < tol) {            
            return it;   
        }
    }
    
    return it;
}

void CSolution::SteadyProblem(Int backend) 
{    
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
    
    // solve for solv.sys.u
    NewtonSolver(solv.sys, disc, prec, backend); 
    
    // update disc.sol.udg
    disc.updateUDG(solv.sys.u, backend); 
}

void CSolution::DIRK(Int backend)
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
        
        //compute time-average solution
        if (disc.common.compudgavg == 1) {
            ArrayAXPBY(disc.sol.udgavg, disc.sol.udgavg, disc.sol.udg, one, one, disc.common.ndofudg1, backend);            
            ArrayAddScalar(&disc.sol.udgavg[disc.common.ndofudg1], one, 1, backend);
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

void CSolution::SolveProblem(Int backend) 
{                  
    this->InitSolution(backend); 
    
    if (disc.common.tdep==1) {        
        // solve time-dependent problems using DIRK
        this->DIRK(backend);            
    }
    else {
        // solve steady-state problems
        this->SteadyProblem(backend);        
        
        // save solutions into binary files            
        this->SaveSolutions(backend);    
        this->SaveSolutionsOnBoundary(backend);         
        if (disc.common.nce>0)
            this->SaveOutputCG(backend);    
    }        
}

//*************** Solver for the PDE model with subproblem **************//

void Residual(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, CSolution *subprob, 
        cublasHandle_t handle, Int backend)
{  
    if (common.subproblem==1) {        
        int npe = common.npe;
        int ne = common.ne;             
        
        // transfer sol.udg to subprob->disc.sol.wdg
        int ncw = subprob->disc.common.ncw;        
        for (int n=0; n<ncw; n++) {            
            int k = subprob->disc.common.vindx[ncw+n]; // solution index
            ArrayExtract(subprob->disc.res.Ru, sol.udg, npe, common.nc, ne, 0, npe, k, k+1, 0, ne, backend);         
            ArrayInsert(subprob->disc.sol.wdg, subprob->disc.res.Ru, npe, ncw, ne, 0, npe, n, n+1, 0, ne, backend);          
        }                
        
        // solving the subproblem to obtain subprob->disc.sol.udg
        subprob->SolveProblem(backend);
        
        // transfer subprob->disc.sol.udg to sol.wdg
        ncw = common.ncw;    
        for (int n=0; n<ncw; n++) {            
            int k = common.vindx[ncw+n]; // solution index
            ArrayExtract(res.Ru, subprob->disc.sol.udg, npe, subprob->disc.common.nc, ne, 0, npe, k, k+1, 0, ne, backend);         
            ArrayInsert(sol.wdg, res.Ru, npe, ncw, ne, 0, npe, n, n+1, 0, ne, backend);          
        }                
    }
        
    Residual(sol, res, app, master, mesh, tmp, common, handle, backend);
}

void Residual(dstype* Ru, dstype* u, CDiscretization& disc, CSolution *subprob, Int backend)
{
    // insert u into udg
    ArrayInsert(disc.sol.udg, u, disc.common.npe, disc.common.nc, disc.common.ne, 0, disc.common.npe, 
            0, disc.common.ncu, 0, disc.common.ne1, backend);  

    // compute the residual vector R(u)
    Residual(disc.sol, disc.res, disc.app, disc.master, disc.mesh, disc.tmp, disc.common, 
            subprob, disc.common.cublasHandle, backend);

    // copy the residual vector to Ru
    ArrayCopy(Ru, disc.res.Ru, disc.common.ndof1, backend);          
}

void MatVec(dstype *w, solstruct &sol, resstruct &res, appstruct &app, masterstruct &master,
      meshstruct &mesh, tempstruct &tmp, commonstruct &common, CSolution *subprob, 
     cublasHandle_t handle, dstype *v, dstype *u, dstype *Ru, Int backend)
{   
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)    
    Int npe = common.npe; // number of nodes on master element    
    Int ne = common.ne1; // number of elements in this subdomain 
    Int nd = common.nd;
    Int N = npe*ncu*ne;
    
    Int order = common.matvecOrder;
    dstype epsilon = common.matvecTol;
#ifdef HAVE_ENZYME
//TODO: there might be a cleaner way to do this...matvecAD v. matvecFD functions? 
    ArrayInsert(sol.dudg, v, npe, nc, ne, 0, npe, 0, ncu, 0, ne, backend); //insert v into dudgs
    dResidual(sol, res, app, master, mesh, tmp, common, handle, backend);
    ArrayAXPBY(w, u, res.dRu, 0.0, 1, N, backend);
#else
    if (order==1) {
        // calculate w = u + epsilon*v
        ArrayAXPBY(w, u, v, 1.0, epsilon, N, backend);
                            
        // insert (u+epsilon*v) into udg
        ArrayInsert(sol.udg, w, npe, nc, ne, 0, npe, 0, ncu, 0, ne, backend);  

        // compute the residual R(u+epsilon*v)
        Residual(sol, res, app, master, mesh, tmp, common, subprob, handle, backend);

        // calculate w = J(u)*v = (R(u+epsilon*v)-R(u))/epsilon    
        ArrayAXPBY(w, res.Ru, Ru, 1.0/epsilon, -1.0/epsilon, N, backend);
    }
    else if (order==2) {
        // calculate w = u - epsilon*v
        ArrayAXPBY(w, u, v, 1.0, -epsilon, N, backend);

        // insert (u-epsilon*v) into udg
        ArrayInsert(sol.udg, w, npe, nc, ne, 0, npe, 0, ncu, 0, ne, backend);  

        // compute the residual R(u-epsilon*v)
        Residual(sol, res, app, master, mesh, tmp, common, subprob, handle, backend);

        // copy res.Ru to Ru
        ArrayCopy(Ru, res.Ru, N, backend);
        
        // calculate w = u + epsilon*v
        ArrayAXPBY(w, u, v, 1.0, epsilon, N, backend);

        // insert (u+epsilon*v) into udg
        ArrayInsert(sol.udg, w, npe, nc, ne, 0, npe, 0, ncu, 0, ne, backend);  

        // compute the residual R(u+epsilon*v)
        Residual(sol, res, app, master, mesh, tmp, common, subprob, handle, backend);
        
        // calculate w = J(u)*v = (R(u+epsilon*v)-R(u-epsilon*v))/(2*epsilon)    
        ArrayAXPBY(w, res.Ru, Ru, 0.5/epsilon, -0.5/epsilon, N, backend);
    }
    else
        error("Matrix-vector multiplication order is not implemented");
#endif
}

void MatVec(dstype *w, dstype *v, dstype *u, dstype *Ru, CDiscretization& disc, CSolution *subprob, Int backend)
{
  MatVec(w, disc.sol, disc.res, disc.app, disc.master, disc.mesh, disc.tmp, disc.common, subprob, 
            disc.common.cublasHandle, v, u, Ru, backend); 
}

void ComputeInitialGuessAndPreconditioner(sysstruct& sys, CDiscretization& disc, 
        CPreconditioner& pre, CSolution *subprob, Int backend)
{           
    Int N = disc.common.ndof1;
    Int RBdim = disc.common.RBcurrentdim;
    dstype *RBcoef = &disc.tmp.tempn[0];
    dstype *RBcoef_tmp = &disc.tmp.tempn[RBdim];
    dstype *H_tmp = &disc.tmp.tempg[RBdim*RBdim];
    
    // compute V = J(u)*W        
    int i = disc.common.RBremovedind -1;
    if (disc.common.RBremovedind==0)
         i = RBdim-1;
    //disc.evalMatVec(&pre.precond.U[i*N], &pre.precond.W[i*N], sys.u, sys.b, backend);      
    MatVec(&pre.precond.U[i*N], &pre.precond.W[i*N], sys.u, sys.b, disc, subprob, backend);
    
    /* RBcoef_tmp = V^T b */
    PGEMTV(disc.common.cublasHandle, N, RBdim, &one, pre.precond.U, N, sys.b, inc1, 
            &zero, RBcoef_tmp, inc1, H_tmp, backend);
    
    /* H = V^T V */
    PGEMTM(disc.common.cublasHandle, RBdim, RBdim, N, &one, pre.precond.U, N, pre.precond.U, N, &zero, 
            pre.precond.H, RBdim, H_tmp, backend);    
    
    /* Solve the linear system */
    Inverse(disc.common.cublasHandle, pre.precond.H, H_tmp, pre.precond.ipiv, RBdim, 1, backend);       

    // RBcoef = inverse(H)*RBcoef_tmp
    PGEMNV(disc.common.cublasHandle, RBdim, RBdim, &one, pre.precond.H, RBdim, 
            RBcoef_tmp, inc1, &zero, RBcoef, inc1, backend);    
        
    // compute x = -W*RBcoef
    PGEMNV(disc.common.cublasHandle, N, RBdim, &minusone, pre.precond.W, N, RBcoef, 
            inc1, &zero, sys.x, inc1, backend);                                                 
}

int GMRES(sysstruct &sys, CDiscretization &disc, CPreconditioner& prec, 
        CSolution *subprob, Int backend)
{
    Int maxit, nrest, orthogMethod, n1, i, k, j = 0;
    Int ncu = disc.common.ncu;
    Int npe = disc.common.npe;
    Int ne = disc.common.ne1;
    Int N = npe*ncu*ne;
    dstype nrmb, nrmr, tol, scalar;
    tol = min(0.01,disc.common.linearSolverTol*disc.common.linearSolverTolFactor);
    maxit = disc.common.linearSolverMaxIter;
    nrest = disc.common.gmresRestart;
    orthogMethod = disc.common.gmresOrthogMethod;
    n1 = nrest + 1;
                
    dstype *s, *y, *cs, *sn, *H;
    s = &sys.tempmem[0];
    y = &sys.tempmem[n1];
    cs = &sys.tempmem[2*n1];
    sn = &sys.tempmem[3*n1];
    H = &sys.tempmem[4*n1];
    
    // calculate Ritz values of polynomial preconditioner
    if  (disc.common.ppdegree > 1) {
        getPoly(disc, prec, sys, sys.lam, sys.randvect, sys.ipiv, N, disc.common.ppdegree, backend);        
        dstype lmin=1.0e10, lmax=-1.0e10, smax=0.0; 
        int m = disc.common.ppdegree;
        for (int i=0; i<m; i++) {
            lmax = (sys.lam[2*m+i] > lmax) ? sys.lam[2*m+i] : lmax;
            lmin = (sys.lam[2*m+i] < lmin) ? sys.lam[2*m+i] : lmin;
            smax = (fabs(sys.lam[3*m+i]) > smax) ? fabs(sys.lam[3*m+i]) : smax;
        }        
        if (smax == 0.0) {
            for (int i=1; i<=m; i++) 
                sys.lam[i-1] = lmin + (lmax-lmin)*(0.5 - 0.5*cos((2*i-1)*M_PI/(2*m))/cos(M_PI/(2*m)));            
            LejaSort(&sys.lam[2*m], &sys.lam[3*m], sys.lam, &sys.lam[3*m], &sys.lam[4*m], m);
        }        
//         for (int i=0; i<disc.common.ppdegree; i++)
//             if (disc.common.mpiRank==0) cout<<sys.lam[2*disc.common.ppdegree+i]<<"  "<<sys.lam[3*disc.common.ppdegree+i]<<endl;
    }    
    
    // compute b
    nrmb = PNORM(disc.common.cublasHandle, N, sys.b, backend);    
    
    scalar = PNORM(disc.common.cublasHandle, N, sys.x, backend);    
    if (scalar>1e-12) {
        // r = A*x
        //disc.evalMatVec(sys.r, sys.x, sys.u, sys.b, backend);
        MatVec(sys.r, sys.x, sys.u, sys.b, disc, subprob, backend);

        // compute the new RHS vector: r = -b - A*x
        ArrayAXPBY(sys.r, sys.b, sys.r, minusone, minusone, N, backend);    

        // norm of the new RHS vector
        nrmr = PNORM(disc.common.cublasHandle, N, sys.r, backend);
    }
    else {
        ArrayAXPBY(sys.r, sys.b, sys.b, minusone, zero, N, backend);
        nrmr = nrmb;  
    }

    if (disc.common.mpiRank==0)
        cout<<"Old RHS Norm: "<<nrmb<<",  New RHS Norm: "<<nrmr<<endl; 
    
    disc.common.linearSolverTolFactor = nrmb/nrmr;
    tol = min(0.1,disc.common.linearSolverTol*disc.common.linearSolverTolFactor);
        
    //disc.common.ppdegree = 0;
    
    // compute r = P*r
    if (disc.common.ppdegree>1) {
        ApplyPreconditioner(sys.r, disc, prec, backend);
        //ApplyComponentNorm(disc, sys.normcu, sys.r, disc.res.Ru, ncu, ncu, npe, ne, backend);
        ApplyPoly(sys.r, disc, prec, sys, sys.q, sys.p, N, backend);    
    }
    else if (disc.common.RBcurrentdim>=0) {
        ApplyPreconditioner(sys.r, disc, prec, backend);
        //ApplyComponentNorm(disc, sys.normcu, sys.r, disc.res.Ru, ncu, ncu, npe, ne, backend);
    }
    
    // compute ||r||
    nrmb = PNORM(disc.common.cublasHandle, N, sys.r, backend);
    nrmr = nrmb;
        
    j = 1;
    while (j < maxit) {
        // v = r/||r||
        ArrayAXPB(sys.v, sys.r, one/nrmr, zero, N, backend);
        
        // initialize s
        s[0] = nrmr;
        for (k=1; k<n1; k++)
            s[k] = zero;
        
        for (i = 0; i < nrest && j < maxit; i++, j++) {
            Int m = i + 1;
            
            // compute v[m] = A*v[i]
            START_TIMING;                        
            //disc.evalMatVec(&sys.v[m*N], &sys.v[i*N], sys.u, sys.b, backend);         
            MatVec(&sys.v[m*N], &sys.v[i*N], sys.u, sys.b, disc, subprob, backend);
            END_TIMING_DISC(0);    
                                    
            START_TIMING;                                    
            // compute v[m] = P*v[m]
            if (disc.common.ppdegree>1) {
                ApplyPreconditioner(&sys.v[m*N], disc, prec, backend);   
                //ApplyComponentNorm(disc, sys.normcu, &sys.v[m*N], disc.res.Ru, ncu, ncu, npe, ne, backend);
                ApplyPoly(&sys.v[m*N], disc, prec, sys, sys.q, sys.p, N, backend);
            }
            else if (disc.common.RBcurrentdim>=0) {
                ApplyPreconditioner(&sys.v[m*N], disc, prec, backend);   
                //ApplyComponentNorm(disc, sys.normcu, &sys.v[m*N], disc.res.Ru, ncu, ncu, npe, ne, backend);
            }
            END_TIMING_DISC(2);    
                        
            // orthogonalize Krylov vectors
            START_TIMING;    
            if (orthogMethod == 0)
                MGS(disc.common.cublasHandle, sys.v, &H[n1*i], N, m, backend);
            else
                CGS(disc.common.cublasHandle, sys.v, &H[n1*i], y, N, m, backend);
            END_TIMING_DISC(1);    
            
            // Apply Givens rotation to compute s
            ApplyGivensRotation(&H[n1*i], s, cs, sn, i, 0);
            
            // compute relative error
            disc.common.linearSolverRelError = fabs(s[i+1])/nrmb;
            
            // check convergence and update solution: x = x + v*s
            if (disc.common.linearSolverRelError < tol) {                
                UpdateSolution(disc.common.cublasHandle, sys.x, y, H, s, sys.v, i, N, n1, backend);
                return j;
            }            
        }        
        
        // update solution: x = x + v*s
        UpdateSolution(disc.common.cublasHandle, sys.x, y, H, s, sys.v, nrest-1, N, n1, backend);
               
        // compute r = A*x
        //disc.evalMatVec(sys.r, sys.x, sys.u, sys.b, backend);
        MatVec(sys.r, sys.x, sys.u, sys.b, disc, subprob, backend);
        
        // r = -b - A*x
        ArrayAXPBY(sys.r, sys.b, sys.r, minusone, minusone, N, backend);
        
        // r = P*r = P*(-b-A*x)
        if (disc.common.ppdegree>1) {
            ApplyPreconditioner(sys.r, disc, prec, backend);            
            //ApplyComponentNorm(disc, sys.normcu, sys.r, disc.res.Ru, ncu, ncu, npe, ne, backend);
            ApplyPoly(sys.r, disc, prec, sys, sys.q, sys.p, N, backend);
        }
        else if (disc.common.RBcurrentdim>=0) {
            ApplyPreconditioner(sys.r, disc, prec, backend);            
            //ApplyComponentNorm(disc, sys.normcu, sys.r, disc.res.Ru, ncu, ncu, npe, ne, backend);
        }
        
        // compute relative error
        nrmr = PNORM(disc.common.cublasHandle, N, sys.r, backend);
        disc.common.linearSolverRelError = nrmr/nrmb;
        
        // check convergence
        if (disc.common.linearSolverRelError < tol) {
            return j;
        }
    }       
    
    if (disc.common.linearSolverRelError > tol) {
        if (disc.common.mpiRank==0) {
            printf("Warning: GMRES(%d) does not converge to the tolerance %g within % d iterations\n",nrest,tol,maxit);
            printf("Warning: The current relative error is %g \n",disc.common.linearSolverRelError);
        }
    }
    return j;
}

int LinearSolver(sysstruct &sys, CDiscretization& disc, CPreconditioner& prec,
         CSolution *subprob, Int backend)
{            
    // evaluate the residual R(u) and set it to sys.b
    //disc.evalResidual(sys.b, sys.u, backend);
    Residual(sys.b, sys.u, disc, subprob, backend);

    int N = disc.common.ndof1;

    // residual norm
    dstype oldnrm = PNORM(disc.common.cublasHandle, N, sys.b, backend); 
                   
    if (oldnrm < disc.common.nonlinearSolverTol) {
      ArraySetValue(sys.x, zero, disc.common.ndof1, backend);
      return 1;           
    }
    
    // construct the preconditioner
    if (disc.common.RBcurrentdim>0) {
        ComputeInitialGuessAndPreconditioner(sys, disc, prec, subprob, backend); 
        
        // v = u + x 
        ArrayAXPBY(sys.v, sys.u, sys.x, one, one, N, backend);  
        //disc.evalResidual(sys.r, sys.v, backend);  
        Residual(sys.r, sys.v, disc, subprob, backend);
        dstype nrmr = PNORM(disc.common.cublasHandle, N, sys.r, backend);
                
        if (nrmr>1.05*oldnrm) {
            ArraySetValue(sys.x, zero, disc.common.ndof1, backend);
            // reset the reduced basis
            disc.common.RBremovedind = 0;
            disc.common.RBcurrentdim = 0;
        }

        if (nrmr < disc.common.nonlinearSolverTol) 
            return 1;               
    }    
    else
        ArraySetValue(sys.x, zero, disc.common.ndof1, backend);
    
    disc.common.linearSolverIter = GMRES(sys, disc, prec, subprob, backend);  
    
    if ((disc.common.linearSolverRelError <= disc.common.linearSolverTol*disc.common.linearSolverTolFactor) && (disc.common.mpiRank==0))             
        printf("GMRES converges to the tolerance %g within % d iterations and %d RB dimensions\n",disc.common.linearSolverTol,disc.common.linearSolverIter,disc.common.RBcurrentdim);        
                 
    return 0;
}

Int NewtonSolver(sysstruct &sys, CDiscretization& disc, CPreconditioner& prec, CSolution *subprob, Int backend)       
{
    Int N = disc.common.ndof1;     
    Int it = 0, maxit = disc.common.nonlinearSolverMaxIter;  
    dstype nrmr, tol;
    tol = disc.common.nonlinearSolverTol; // tolerance for the residual
        
    nrmr = PNORM(disc.common.cublasHandle, N, sys.u, backend);
    if (disc.common.mpiRank==0)
        cout<<"Newton Iteration: "<<it<<",  Solution Norm: "<<nrmr<<endl;                                                    
    
    // use PTC to solve the system: R(u) = 0
    for (it=1; it<maxit; it++) {                        
        //nrmrold = nrmr;
        
        // solve the linear system: (lambda*B + J(u))x = -R(u)
        int status;
        status = LinearSolver(sys, disc, prec, subprob, backend);
                                
        // update the solution: u = u + x
        ArrayAXPBY(sys.u, sys.u, sys.x, one, one, N, backend); 

        if (status==1)
           return it;
       
        // compute both the residual vector and sol.udg  
        //disc.evalResidual(sys.r, sys.u, backend);
        Residual(sys.r, sys.u, disc, subprob, backend);
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
        }
                
        if (disc.common.mpiRank==0)
            cout<<"Newton Iteration: "<<it<<",  Residual Norm: "<<nrmr<<endl;                           
                        
        // update the reduced basis
        if (status==0)
            UpdateRB(sys, disc, prec, backend);      
        
        // check convergence
        if (nrmr < tol) {            
            return it;   
        }
    }
    
    return it;
}

void CSolution::SteadyProblem(CSolution *subprob, Int backend) 
{    
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
    
    // solve for solv.sys.u
    NewtonSolver(solv.sys, disc, prec, subprob, backend); 
    
    // update disc.sol.udg
    disc.updateUDG(solv.sys.u, backend); 
    
    // solve the subproblem
    if (disc.common.subproblem==1) {        
        int npe = disc.common.npe;
        int ne = disc.common.ne;             
        
        // transfer sol.udg to subprob->disc.sol.wdg
        int ncw = subprob->disc.common.ncw;        
        for (int n=0; n<ncw; n++) {            
            int k = subprob->disc.common.vindx[ncw+n]; // solution index
            ArrayExtract(subprob->disc.res.Ru, disc.sol.udg, npe, disc.common.nc, ne, 0, npe, k, k+1, 0, ne, backend);         
            ArrayInsert(subprob->disc.sol.wdg, subprob->disc.res.Ru, npe, ncw, ne, 0, npe, n, n+1, 0, ne, backend);          
        }                
        
        // solving the subproblem to obtain subprob->disc.sol.udg
        subprob->SolveProblem(backend);
        
        // transfer subprob->disc.sol.udg to sol.wdg
        ncw = disc.common.ncw;        
        for (int n=0; n<ncw; n++) {            
            int k = disc.common.vindx[ncw+n]; // solution index
            ArrayExtract(disc.res.Ru, subprob->disc.sol.udg, npe, subprob->disc.common.nc, ne, 0, npe, k, k+1, 0, ne, backend);         
            ArrayInsert(disc.sol.wdg, disc.res.Ru, npe, ncw, ne, 0, npe, n, n+1, 0, ne, backend);          
        }                
    }    
}

void CSolution::DIRK(CSolution *subprob, Int backend)
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
            this->SteadyProblem(subprob, backend);                             

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
            ArrayAXPBY(disc.sol.udgavg, disc.sol.udgavg, disc.sol.udg, one, one, disc.common.ndofudg1, backend);            
            ArrayAddScalar(&disc.sol.udgavg[disc.common.ndofudg1], one, 1, backend);
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

void CSolution::SolveProblem(CSolution *subprob, Int backend) 
{                  
    this->InitSolution(backend); 
            
    if (disc.common.tdep==1) {        
        // solve time-dependent problems using DIRK
        this->DIRK(subprob, backend);            
    }
    else {
        // solve steady-state problems
        this->SteadyProblem(subprob, backend);        
        
        // save solutions into binary files            
        this->SaveSolutions(backend);    
        this->SaveSolutionsOnBoundary(backend);         
        if (disc.common.nce>0)
            this->SaveOutputCG(backend);    
    }        
}
