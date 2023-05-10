#ifndef __SOLUTION
#define __SOLUTION

#include "solution.h"
#include "previoussolutions.cpp"
#include "savesolutions.cpp"
#include "updatesolution.cpp"
#include "updatesource.cpp"
#include "timestepcoeff.cpp"
#include "subproblem.cpp"

#ifdef TIMESTEP  
#include <sys/time.h>
#endif

/// Q: why DG2CG2 and all of the communication twice?
////A: sometimes needed, write a loop
/// Q: there are two seperate GetQ and GetW, for interior and edges...where should this go in between here?
///   after all getQ's or is there a more intelligent way it can be done? 
/// A: do it in parts for interior and interface
/// Q: communication is actually slightly different; we send (u, q) while in residual appears to be only (u)
///    should we just send u and calculate q as needed? or resend information...
/// A: dont need to send q, q is recalculated. DO need to send/receive odg
/// Q: can avfield and odg be used at the same time? Can ODG be put in to residual evaluation? 
/// A: not now
/// Q: about general MPI residual: need to send udg and dudg; repeat communication or append on to end of buffsend? 
/// A: load in same buffsend, if careful about it

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
            ArrayExtract(disc.res.Ru, disc.sol.udgarray[m], npe, disc.common.ncarray[m], ne, 0, npe, k, k+1, 0, ne, backend);         
            // insert Ru into odg
            ArrayInsert(disc.sol.odg, disc.res.Ru, npe, nco, ne, 0, npe, n, n+1, 0, ne, backend);          
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
            GetArrayAtIndex(disc.tmp.buffsend, disc.sol.odg, disc.mesh.elemsendodg, bsz*disc.common.nelemsend, backend);

#ifdef HAVE_CUDA
            cudaDeviceSynchronize();
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
            PutArrayAtIndex(disc.sol.odg, disc.tmp.buffrecv, disc.mesh.elemrecvodg, bsz*disc.common.nelemrecv, backend);
#endif
    //    END_TIMING_DISC(98);    
        }

    }
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
    
    // START_TIMING;
    // use PTC to solve steady problem
    solv.PseudoTransientContinuation(disc, prec, out, backend);           
    // END_TIMING_DISC(97);    
        
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
    // compute the geometry quantities
    disc.compGeometry(backend);
    
    //printArray2D(&disc.sol.elemg[0],disc.common.nge,10,backend);
    
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
    // copy sol.dudg to sys.v; should get overwritten if matvec is used by solver   
    ArrayExtract(solv.sys.v, disc.sol.dudg, disc.common.npe, disc.common.nc, disc.common.ne1, 
        0, disc.common.npe, 0, disc.common.ncu, 0, disc.common.ne1, backend);                                               
        
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
                ArraySetValue(disc.sol.udgavg, zero, disc.common.ndofudg1+1, backend);
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

void CSolution::SolveProblemRJv(ofstream &out, Int backend) 
{          
    this->InitSolution(backend); 

    // calls R(U) with mesh.udg, stored in solv.sys.b
    disc.evalResidual(solv.sys.b, solv.sys.u, backend);

    // calls dR/dU(U) v for v loaded in solv.sys.v
    disc.evalMatVec(solv.sys.v, solv.sys.v, solv.sys.u, solv.sys.b, backend); 

    string filenameR = disc.common.fileout + "_Ru_test_np" + NumberToString(disc.common.mpiRank) + ".bin";        
    string filenameJv = disc.common.fileout + "_Jv_test_np" + NumberToString(disc.common.mpiRank) + ".bin";                    
           
    writearray2file(filenameR, solv.sys.b, disc.common.ndof1, backend);
    writearray2file(filenameJv, solv.sys.v, disc.common.ndof1, backend);
}

void CSolution::SaveSolutions(Int backend) 
{
   if (disc.common.tdep==1) { 
        if (((disc.common.currentstep+1) % disc.common.saveSolFreq) == 0)             
        {        
            string filename = disc.common.fileout + "_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank) + ".bin";     
            if (disc.common.saveSolOpt==0)
                writearray2file(filename, solv.sys.u, disc.common.ndof1, backend);
            else
                writearray2file(filename, disc.sol.udg, disc.common.ndofudg1, backend);
                        
            if (disc.common.ncw>0) {
                string fn = disc.common.fileout + "_wdg_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank) + ".bin";                    
                writearray2file(fn, solv.sys.wtmp, disc.common.ndofw1, backend);
            }                        
            
            if (disc.common.compudgavg == 1) {
                string fn1 = disc.common.fileout + "avg_np" + NumberToString(disc.common.mpiRank) + ".bin"; 
                writearray2file(fn1, disc.sol.udgavg, disc.common.ndofudg1+1, backend);
            }
        }    
   }
   else {
        string filename = disc.common.fileout + "_np" + NumberToString(disc.common.mpiRank) + ".bin";                    
        if (disc.common.saveSolOpt==0)
            writearray2file(filename, solv.sys.u, disc.common.ndof1, backend);
        else
            writearray2file(filename, disc.sol.udg, disc.common.ndofudg1, backend);       
        
   }    
}

void CSolution::ReadSolutions(Int backend) 
{
   if (disc.common.tdep==1) { 
        if (((disc.common.currentstep+1) % disc.common.saveSolFreq) == 0)             
        {        
            string filename = disc.common.fileout + "_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank) + ".bin";     
            if (disc.common.saveSolOpt==0) {
                readarrayfromfile(filename, &solv.sys.u, disc.common.ndof1, backend);
                // insert u into udg
                ArrayInsert(disc.sol.udg, solv.sys.u, disc.common.npe, disc.common.nc, 
                 disc.common.ne, 0, disc.common.npe, 0, disc.common.ncu, 0, disc.common.ne1, backend);  
            }
            else
                readarrayfromfile(filename, &disc.sol.udg, disc.common.ndofudg1, backend);        
            
            if (disc.common.ncw>0) {
                string fn = disc.common.fileout+"_wdg_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank) + ".bin";                    
                readarrayfromfile(fn, &disc.sol.wdg, disc.common.ndofw1, backend);     
            }                                    
        }                                
   }
   else {
        string filename = disc.common.fileout + "_np" + NumberToString(disc.common.mpiRank) + ".bin";                    
        if (disc.common.saveSolOpt==0) {
            readarrayfromfile(filename, &solv.sys.u, disc.common.ndof1, backend);
            // insert u into udg
            ArrayInsert(disc.sol.udg, solv.sys.u, disc.common.npe, disc.common.nc, 
             disc.common.ne, 0, disc.common.npe, 0, disc.common.ncu, 0, disc.common.ne1, backend);              
        }
        else
            readarrayfromfile(filename, &disc.sol.udg, disc.common.ndofudg1, backend);      
        
        if (disc.common.ncw>0) {
            string fn = disc.common.fileout + "_wdg_np" + NumberToString(disc.common.mpiRank) + ".bin";                    
            readarrayfromfile(fn, &disc.sol.wdg, disc.common.ndofw1, backend);     
        }                                            
   }    
}
 
void CSolution::SaveOutputDG(Int backend) 
{
   if (disc.common.tdep==1) { 
        if (((disc.common.currentstep+1) % disc.common.saveSolFreq) == 0)             
        {                    
            string filename1 = disc.common.fileout + "_outputDG_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank) + ".bin";     
            disc.evalOutput(solv.sys.v, backend);                        
            writearray2file(filename1, solv.sys.v, disc.common.ndofedg1, backend);       
        }                                
   }
   else {
        string filename1 = disc.common.fileout + "_outputDG_np" + NumberToString(disc.common.mpiRank) + ".bin";                           
        disc.evalOutput(solv.sys.v, backend);
        writearray2file(filename1, solv.sys.v, disc.common.ndofedg1, backend);       
   }    
}

void CSolution::SaveOutputCG(Int backend) 
{
   if (disc.common.tdep==1) { 
        if (((disc.common.currentstep+1) % disc.common.saveSolFreq) == 0)             
        {                    
            string filename1 = disc.common.fileout + "_outputCG_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank) + ".bin";     
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
        string filename1 = disc.common.fileout + "_outputCG_np" + NumberToString(disc.common.mpiRank) + ".bin";                            
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
                    GetArrayAtIndex(disc.tmp.tempn, disc.sol.udg, &disc.mesh.findudg1[npf*nc*f1], nn*nc, backend);
                    string filename = disc.common.fileout + "bou_t" + NumberToString(disc.common.currentstep+disc.common.timestepOffset+1) + "_np" + NumberToString(disc.common.mpiRank) + ".bin";     
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
                GetArrayAtIndex(disc.tmp.tempn, disc.sol.xdg, &disc.mesh.findxdg1[npf*ncx*f1], nn*ncx, backend);
                string filename = disc.common.fileout + "node_np" + NumberToString(disc.common.mpiRank) + ".bin";     
                writearray2file(filename, disc.tmp.tempn, nn*ncx, backend);            
            }
        }                                                                   
    }
}

#endif        
