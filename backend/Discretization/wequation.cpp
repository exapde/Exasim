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

#endif

