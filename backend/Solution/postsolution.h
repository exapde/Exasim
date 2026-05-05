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
#ifndef __POSTSOLUTION_H__
#define __POSTSOLUTION_H__

#include <exasim/detail/abi_adapter.hpp>

template <class M = exasim::detail::AbiAdapter>
class CSolution {
private:
public:
    CDiscretization<M> disc;  // spatial discretization class
    CPreconditioner<M> prec;  // precondtioner class
    CSolver<M> solv;          // linear and nonlinear solvers
    CVisualization vis;    // visualization class
    std::ofstream outsol;       // storing solutions
    std::ofstream outwdg;  
    std::ofstream outuhat;  
    std::ofstream outbouxdg;  
    std::ofstream outboundg;  
    std::ofstream outbouudg;  
    std::ofstream outbouwdg;  
    std::ofstream outbouuhat;    
    std::ofstream outqoi;
    
    // constructor 
    CSolution(std::string filein, std::string fileout, std::string exasimpath, Int mpiprocs, Int mpirank, 
              Int fileoffset, Int omprank, Int backend, Int builtinmodelID,
              Int nsca, Int nvec, Int nten, Int nsurf, Int nvqoi)   
       : disc(filein, fileout, exasimpath, mpiprocs, mpirank, fileoffset, omprank, backend, builtinmodelID, nsca, nvec, nten, nsurf, nvqoi),
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

        if ((disc.common.nintfaces > 0) && (disc.common.coupledcondition>0)) disc.common.ne0 = disc.common.intepartpts[0];            

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
    };        
    
    // destructor        
    ~CSolution() { 
        if (outqoi.is_open()) { outqoi.close(); }
    }; 

    // precompute some quantities
    void InitSolution(Int backend);    
        
    // save solutions in binary files
    void SaveSolutions(Int backend);    

    void SaveQoI(Int backend);    

    // save fields in VTU files
    void SaveParaview(Int backend, std::string fname_modifier = "", bool force_tdep_write = false);    
    
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


