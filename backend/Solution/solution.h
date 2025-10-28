/**
 * @file solution.h
 * @brief Defines the CSolution class for managing and solving PDE problems.
 *
 * The CSolution class encapsulates the main components and routines for
 * discretization, preconditioning, and solving linear/nonlinear systems
 * arising from PDEs. It provides interfaces for steady-state and time-dependent
 * problem solving, solution initialization, and input/output operations.
 */

 /**
  * @class CSolution
  * @brief Main class for handling PDE solution workflow.
  *
  * This class manages the discretization, preconditioning, and solver modules,
  * and provides methods for solving steady-state and time-dependent problems,
  * as well as saving and loading solutions and outputs.
  *
  * @section Members
  * - CDiscretization disc: Handles spatial discretization.
  * - CPreconditioner prec: Manages preconditioning for solvers.
  * - CSolver solv: Provides linear and nonlinear solver routines.
  *
  * @section Methods
  * - CSolution(...): Constructor initializing all components.
  * - ~CSolution(): Destructor.
  * - SteadyProblem(...): Solve steady-state problems.
  * - SteadyProblem_PTC(...): Solve steady-state problems using PTC.
  * - TimeStepping(...): Advance solution in time using DIRK/BDF schemes.
  * - UnsteadyProblem(...): Solve time-dependent problems.
  * - DIRK(...): Time integration using DIRK scheme.
  * - InitSolution(...): Precompute quantities for solution.
  * - SolveProblem(...): Solve both steady-state and time-dependent problems.
  * - SaveSolutions(...): Save solutions to binary files.
  * - SaveSolutionsOnBoundary(...): Save boundary solutions to binary files.
  * - SaveNodesOnBoundary(...): Save boundary nodes to binary files.
  * - ReadSolutions(...): Read solutions from binary files.
  * - SaveOutputDG(...): Save DG output to binary files.
  * - SaveOutputCG(...): Save CG output to binary files.
  */
#ifndef __SOLUTION_H__
#define __SOLUTION_H__

// Common helper: open file and write 3-element header [a0, a1, a2]
void open_and_write(std::ofstream& ofs,
                    const std::string& prefix,
                    int rank, int offset,
                    int a0, int a1, int a2,
                    const std::string& fileout)
{
    std::string filename = fileout + prefix +
                           NumberToString(rank - offset) + ".bin";
    ofs.open(filename.c_str(), std::ios::out | std::ios::binary);
    if (!ofs) error("Failed to open file: " + filename);

    dstype a[3] = { dstype(a0), dstype(a1), dstype(a2) };
    writearray(ofs, a, 3);
}

class CSolution {
private:
public:
    CDiscretization disc;  // spatial discretization class
    CPreconditioner prec;  // precondtioner class 
    CSolver solv;          // linear and nonlinear solvers
    CVisualization vis;    // visualization class
    ofstream outsol;       // storing solutions
    ofstream outwdg;  
    ofstream outuhat;  
    ofstream outbouxdg;  
    ofstream outboundg;  
    ofstream outbouudg;  
    ofstream outbouwdg;  
    ofstream outbouuhat;  
    ofstream outqoi;
    
    // constructor 
    CSolution(string filein, string fileout, Int mpiprocs, Int mpirank, Int fileoffset, Int omprank, Int backend)   
       : disc(filein, fileout, mpiprocs, mpirank, fileoffset, omprank, backend),
         prec(disc, backend), solv(disc, backend), vis(disc, backend) 
    {   
        int ncx = disc.common.ncx;                            
        int nd = disc.common.nd;     
        int ncu = disc.common.ncu;     
        int nc = (disc.common.saveSolOpt==0) ? disc.common.ncu : disc.common.nc;        
        int ncw = disc.common.ncw;
        int npe = disc.common.npe;
        int npf = disc.common.npf;
        int ne = disc.common.ne1;     
        int nf = disc.common.nf;     
        int rank = disc.common.mpiRank;
        int offset = disc.common.fileoffset;
        std::string base = disc.common.fileout;

        if (disc.common.spatialScheme > 0) {     
            TemplateFree(disc.res.Minv2, disc.common.backend);
            disc.res.Minv2 = nullptr;
            disc.res.szMinv2 = 0;

            TemplateFree(disc.res.Mass2, disc.common.backend);     
            disc.res.Mass2 = nullptr;
            disc.res.szMass2 = 0;

            if (disc.common.nsurf == 0) {
                TemplateFree(disc.sol.faceg, disc.common.backend);                    
                disc.sol.faceg = nullptr;
                disc.sol.szfaceg = 0;
            }

            if  (disc.common.ppdegree <= 1) {
                TemplateFree(solv.sys.randvect, disc.common.backend);                    
                solv.sys.randvect = nullptr;
                solv.sys.szrandvect = 0;
            }
        }

        if (disc.common.mpiRank==0) {     
            int n1 = disc.app.sizeoffloat();
            int n2 = disc.master.sizeoffloat();
            int n4 = disc.tmp.sizeoffloat();  
            int n5 = disc.res.sizeoffloat();        
            int n6 = disc.sol.sizeoffloat();        
            int n7 = solv.sys.sizeoffloat();   
            int n8 = prec.precond.sizeoffloat();   
            cout<<"Number of doubles for each data structure on rank 0\n";      
            cout<<"app     master     tmp    res    sol    sys     prec\n";             
            printf("%d   %d    %d    %d    %d     %d    %d\n", n1, n2, n4, n5, n6, n7, n8);
            disc.tmp.printinfo();
            disc.res.printinfo();
            disc.sol.printinfo();
            solv.sys.printinfo();
            prec.precond.printinfo();                  
        }

        #ifdef HAVE_HIP
            checkHipMemory(disc.common.mpiRank);  
        #endif
        #ifdef HAVE_CUDA
            checkCudaMemory(disc.common.mpiRank);  
        #endif          

        if (mpirank==0 && (disc.common.nvqoi > 0 || disc.common.nsurf > 0)) {
            outqoi.open(base + "qoi.txt", std::ios::out);                         
            outqoi << std::setw(16) << std::left << "Time";
            for (size_t i = 0; i < disc.common.nvqoi; ++i) {                
                outqoi << std::setw(16) << std::left << "Domain_QoI" + std::to_string(i + 1);
            }
            for (size_t i = 0; i < disc.common.nsurf; ++i) {                
                outqoi << std::setw(16) << std::left << "Boundary_QoI" + std::to_string(i + 1);
            }
            outqoi << "\n";
        }

        open_and_write(outsol, "udg_np", rank, offset, npe, nc, ne, base);

        if (ncw > 0) {     
            open_and_write(outwdg, "wdg_np", rank, offset, npe, ncw, ne, base);
        }

        if (disc.common.spatialScheme==1) {         
            open_and_write(outuhat, "uhat_np", rank, offset, ncu, npf, nf, base);
        }

        if ( disc.common.saveSolBouFreq>0 ) {
            Int nfbou = 0;
            for (Int j=0; j<disc.common.nbf; j++) {
                Int f1 = disc.common.fblks[3*j]-1;
                Int f2 = disc.common.fblks[3*j+1];    
                Int ib = disc.common.fblks[3*j+2];            
                if (ib == disc.common.ibs) {     
                    Int nf = f2-f1;
                    nfbou += nf;
                }
            }

            open_and_write(outbouxdg, "bouxdg_np", rank, offset, npf, nfbou, ncx, base);
            open_and_write(outboundg, "boundg_np", rank, offset, npf, nfbou, nd, base);
            open_and_write(outbouudg, "bouudg_np", rank, offset, npf, nfbou, disc.common.nc, base);            
            open_and_write(outbouuhat, "bouuhat_np", rank, offset, npf, nfbou, ncu, base);
            if (ncw > 0) open_and_write(outbouwdg, "bouwdg_np", rank, offset, npf, nfbou, ncw, base);
        }
        
        // if (!outsol) error("Unable to open file " + filename);        
//           disc.common.printinfo();
//           disc.app.printinfo();
//           disc.res.printinfo();
//           disc.tmp.printinfo();
//           disc.sol.printinfo();
//           disc.mesh.printinfo();
//           disc.master.printinfo();       
//           prec.precond.printinfo();
//           solv.sys.printinfo();
//           error("here");
    };        
    
    // destructor        
    ~CSolution() { 
        if (outsol.is_open()) { outsol.close(); }
        if (outwdg.is_open()) { outwdg.close(); }
        if (outuhat.is_open()) { outuhat.close(); }
        if (outbouxdg.is_open()) { outbouxdg.close(); }
        if (outboundg.is_open()) { outboundg.close(); }
        if (outbouudg.is_open()) { outbouudg.close(); }
        if (outbouwdg.is_open()) { outbouwdg.close(); }
        if (outbouuhat.is_open()) { outbouuhat.close(); }
        if (outqoi.is_open()) { outqoi.close(); }
    }; 

    // solve steady-state problems
    void SteadyProblem(Int backend);
    
    void SteadyProblem(ofstream &out, Int backend);    

    void SteadyProblem_PTC(ofstream &out, Int backend);    
            
    void SteadyProblem(CSolution *subprob, Int backend);    
    
    // advance the solution to the next time step using DIRK and BDF schemes
    void TimeStepping(ofstream &out, dstype time, Int istep, Int backend);
        
    // solve time-dependent problems
    void UnsteadyProblem(ofstream &out, Int backend);    
        
    // solve time-dependent problems
    void DIRK(Int backend);
    
    void DIRK(ofstream &out, Int backend);    
    
    void DIRK(CSolution *subprob, Int backend);

    // precompute some quantities
    void InitSolution(Int backend);    
    
    // solve both steady-state and time-dependent problems
    void SolveProblem(Int backend);    
    
    void SolveProblem(ofstream &out, Int backend);    
    
    void SolveProblem(CSolution *subprob, Int backend);    
    
    // save solutions in binary files
    void SaveSolutions(Int backend);    

    void SaveQoI(Int backend);    

    // save fields in VTU files
    void SaveParaview(Int backend);    
    
    // save solutions in binary files
    void SaveSolutionsOnBoundary(Int backend);    
    
    // save solutions in binary files
    void SaveNodesOnBoundary(Int backend);    
    
    // read solutions from binary files
    void ReadSolutions(Int backend);        
    
    // save output in binary files
    void SaveOutputDG(Int backend);        
    
    // save output in binary files
    void SaveOutputCG(Int backend);        
};

#endif        


