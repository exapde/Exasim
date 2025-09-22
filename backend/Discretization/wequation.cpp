/*
  wequation.cpp

  This file contains functions for solving the auxiliary variable equation (w-equation) in the context of a discontinuous Galerkin (DG) or hybridized DG (HDG) discretization. The w-equation is typically used for time-dependent or differential-algebraic equation (DAE) systems, where auxiliary variables are introduced to facilitate the solution process.

  Functions:

  1. void wEquation(dstype *wdg, dstype *xdg, dstype *udg, dstype *odg, dstype *wsrc, dstype *tempg, appstruct &app, commonstruct &common, Int ng, Int backend)
    - Solves for the auxiliary variable w (wdg) given the solution variables (udg), coordinates (xdg), other variables (odg), and source term (wsrc).
    - Handles both wave-type equations (explicit update) and general DAE systems (Newton iteration).
    - Uses temporary storage (tempg) for intermediate calculations.
    - Supports up to five auxiliary variables (ncw <= 5).

  2. void wEquation(dstype *wdg, dstype *wdg_udg, dstype *xdg, dstype *udg, dstype *odg, dstype *wsrc, dstype *tempg, appstruct &app, commonstruct &common, Int ng, Int backend)
    - Extended version of wEquation that also computes the sensitivity of w with respect to udg (wdg_udg).
    - Uses Newton iteration for nonlinear systems and computes the Jacobian of the source term.
    - Supports up to five auxiliary variables (ncw <= 5).

  3. void GetW(dstype *w, solstruct &sol, tempstruct &tmp, appstruct &app, commonstruct &common, Int backend)
    - Loops over element blocks and computes the auxiliary variable w for each element.
    - Extracts element-wise data from global arrays, solves the w-equation, and writes the results back to the global array.

  Notes:
  - The functions rely on several utility routines for matrix operations, source term evaluation, and error handling.
  - The Newton iteration in wEquation checks for convergence and throws an error if the solution does not converge within the specified tolerance.
  - The code is designed for flexibility in the number of auxiliary variables and supports both explicit and implicit time integration schemes.
*/
#ifndef __WEQUATION
#define __WEQUATION

void wEquation(dstype *wdg, dstype *xdg, dstype *udg, dstype *odg, dstype *wsrc, 
      dstype *tempg, appstruct &app, commonstruct &common, Int ng, Int backend)
{        
    Int ncu = common.ncu; // number of compoments of (u)
    Int nd = common.nd; // spatial dimension
    Int nc = common.nc; // number of compoments of (u, q)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    //Int npe = common.npe; // number of nodes on master element    
    Int modelnumber = common.modelnumber;
    dstype time = common.time;                
    dstype *uinf = app.uinf;
    dstype *physicsparam = app.physicsparam;

    if (common.wave==1) {
        // dw/dt = u -> (dtfactor * w - wsrc) = u -> w = (1/dtfactor) * (u + wsrc)
        dstype scalar = one/common.dtfactor;
        ArrayAXPBY(wdg, udg, wsrc, scalar, scalar, ng*ncw);
    }        
    else {         
        // fix bug here
        dstype *s = tempg; // temporary array    
        dstype *s_wdg = &tempg[ng*ncw]; // temporary array 
        // use Newton to solve the nonlinear system  alpha * dw/dt + beta w = S(w, u, q) to obtain w for given (u, q)             
        for (int iter=0; iter<20; iter++) {
          // alpha * dw/dt + beta w = sourcew(u,q,w) -> alpha (dtfactor * w - wsrc) + beta w = sourcew(u,q,w) 
          // ->  (alpha * dtfactor + beta) w - alpha * wsrc - sourcew(u,q,w) = 0              

          // calculate the source term Sourcew(xdg, udg, odg, wdg)
          HdgSourcewonly(s, s_wdg, xdg, udg, odg, wdg, uinf, physicsparam, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);            
          
          // alpha*dirkd/dt + beta
          dstype scalar = common.dae_alpha*common.dtfactor + common.dae_beta;

          // calculate residual vector = sourcew(u,q,w) + alpha * wsrc - (alpha * dtfactor + beta) w
          ArrayAdd3Vectors(s, s, wsrc, wdg, one, common.dae_alpha, -scalar, ng*ncw);                           

          // compute jacobian matrix = (alpha * dtfactor + beta) - s_wdg
          ArrayAXPB(s_wdg, s_wdg, minusone, scalar, ng*ncw*ncw);                

          // solve the linear system jacobian * dw = residual
          if (ncw==1) {                
            SmallMatrixSolve11(s, s_wdg, ng, 1);
          }
          else if (ncw==2) {
            SmallMatrixSolve22(s, s_wdg, ng, 1);
          }
          else if (ncw==3) {
            SmallMatrixSolve33(s, s_wdg, ng, 1);
          }
          else if (ncw==4) {
            SmallMatrixSolve44(s, s_wdg, ng, 1);
          }
          else if (ncw==5) {
            SmallMatrixSolve55(s, s_wdg, ng, 1);
          }
          else {
            error("DAE functionality supports at most five variables.");
          }              

          // update w = w + dw
          ArrayAXPBY(wdg, wdg, s, one, one, ng*ncw);

          // check convergence
          dstype nrm = NORM(common.cublasHandle, ng*ncw, s, backend);
          if (nrm < 1e-6) break;              
        }                        
    }            
}

void wEquation(dstype *wdg, dstype *wdg_udg, dstype *xdg, dstype *udg, dstype *odg, dstype *wsrc, 
       dstype *tempg, appstruct &app, commonstruct &common, Int ng, Int backend)
{        
    Int ncu = common.ncu; // number of compoments of (u)
    Int nd = common.nd; // spatial dimension
    Int nc = common.nc; // number of compoments of (u, q)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    //Int npe = common.npe; // number of nodes on master element    
    Int modelnumber = common.modelnumber;
    dstype time = common.time;                
    dstype *uinf = app.uinf;
    dstype *physicsparam = app.physicsparam;

    if (common.wave==1) {
        // dw/dt = u -> (dtfactor * w - wsrc) = u -> w = (1/dtfactor) * (u + wsrc)
        dstype scalar = one/common.dtfactor;
        ArrayAXPBY(wdg, udg, wsrc, scalar, scalar, ng*ncw);
        ArraySetValue(wdg_udg, scalar, ng*ncw*nc);
    }        
    else {   
         // fix bug here
        dstype *s = tempg; // temporary array    
        dstype *s_wdg = &tempg[ng*ncw]; // temporary array              
        // use Newton to solve the nonlinear system  alpha * dw/dt + beta w = S(w, u, q) to obtain w for given (u, q)                
        for (int iter=0; iter<20; iter++) {
          // alpha * dw/dt + beta w = sourcew(u,q,w) -> alpha (dtfactor * w - wsrc) + beta w = sourcew(u,q,w) 
          // ->  (alpha * dtfactor + beta) w - alpha * wsrc - sourcew(u,q,w) = 0              

          // calculate the source term Sourcew(xdg, udg, odg, wdg)
          HdgSourcewonly(s, s_wdg, xdg, udg, odg, wdg, uinf, physicsparam, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);            
          
          // alpha*dirkd/dt + beta
          dstype scalar = common.dae_alpha*common.dtfactor + common.dae_beta;

          // calculate residual vector = sourcew(u,q,w) + alpha * wsrc - (alpha * dtfactor + beta) w
          ArrayAdd3Vectors(s, s, wsrc, wdg, one, common.dae_alpha, -scalar, ng*ncw);                           

          // compute jacobian matrix = (alpha * dtfactor + beta) - s_wdg
          ArrayAXPB(s_wdg, s_wdg, minusone, scalar, ng*ncw*ncw);                

          // solve the linear system jacobian * dw = residual
          if (ncw==1) {                
            SmallMatrixSolve11(s, s_wdg, ng, 1);
          }
          else if (ncw==2) {
            SmallMatrixSolve22(s, s_wdg, ng, 1);
          }
          else if (ncw==3) {
            SmallMatrixSolve33(s, s_wdg, ng, 1);
          }
          else if (ncw==4) {
            SmallMatrixSolve44(s, s_wdg, ng, 1);
          }
          else if (ncw==5) {
            SmallMatrixSolve55(s, s_wdg, ng, 1);
          }
          else {
            error("DAE functionality supports at most five variables.");
          }              

          // update w = w + dw
          ArrayAXPBY(wdg, wdg, s, one, one, ng*ncw);

          // check convergence
          dstype nrm = NORM(common.cublasHandle, ng*ncw, s, backend);
          if (nrm < 1e-8) {     
            // wdg_udg is actually s_udg 
            HdgSourcew(s, wdg_udg, s_wdg, xdg, udg, odg, wdg, uinf, physicsparam, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);            
            
            // fix bug here
            // compute jacobian matrix = (alpha * dtfactor + beta) - s_wdg
            ArrayAXPB(s_wdg, s_wdg, minusone, scalar, ng*ncw*ncw);                
            
            // w_udg = -inverse(s_wdg) * s_udg
            if (ncw==1) {                
              SmallMatrixSolve11(wdg_udg, s_wdg, ng, nc);
            }
            else if (ncw==2) {
              SmallMatrixSolve22(wdg_udg, s_wdg, ng, nc);
            }
            else if (ncw==3) {
              SmallMatrixSolve33(wdg_udg, s_wdg, ng, nc);
            }
            else if (ncw==4) {
              SmallMatrixSolve44(wdg_udg, s_wdg, ng, nc);
            }
            else if (ncw==5) {
              SmallMatrixSolve55(wdg_udg, s_wdg, ng, nc);
            }
            else {
              error("DAE functionality supports at most three dependent variables.");
            }                          
            break;              
          }
          else {
            if (iter==20) {
              error("Newton in wequation does not converge to 1e-8.");
            }
          } 
        }                        
    }            
}

void GetW(dstype *w, solstruct &sol, tempstruct &tmp, appstruct &app, commonstruct &common, Int backend)
{
  for (Int j=0; j<common.nbe; j++) {         
      Int e1 = common.eblks[3*j]-1;
      Int e2 = common.eblks[3*j+1];
      Int ns = e2-e1;        
      Int ng = common.npe*ns;
      Int ncw = common.ncw;
      Int ncx = common.ncx;
      Int nc = common.nc;
      Int nco = common.nco;
      dstype* wdg = &tmp.tempn[0];
      dstype* xdg = &tmp.tempn[ng*ncw];
      dstype* udg = &tmp.tempn[ng*(ncw+ncx)];
      dstype* odg = &tmp.tempn[ng*(ncw+ncx+nc)];
      dstype* sdg = &tmp.tempn[ng*(ncw+ncx+nc+nco)];
      GetElemNodes(wdg, w, common.npe, ncw, 0, ncw, e1, e2);
      GetElemNodes(xdg, sol.xdg, common.npe, ncx, 0, ncx, e1, e2);
      GetElemNodes(udg, sol.udg, common.npe, nc, 0, nc, e1, e2);
      GetElemNodes(odg, sol.odg, common.npe, nco, 0, nco, e1, e2);
      GetElemNodes(sdg, sol.wsrc, common.npe, ncw, 0, ncw, e1, e2);
      wEquation(wdg, xdg, udg, odg, sdg, tmp.tempg, app, common, ng, common.backend);
      PutElemNodes(w, wdg, common.npe, ncw, 0, ncw, e1, e2);
  }   
}

#endif

