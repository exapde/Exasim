function [Un,Hn,Pn] = hdgdirk(master,mesh,pde,UDG,UH,PDG,time,dt,nstage,torder)
%HDG_SOLVE_DIRK Solve using the HDG method and Newton's iteraion - DIRK timestepping
%   [UH,QH,UHAT] = HDG_SOLVE_DIRK(MASTER,MESH,UH,QH,UHAT,APP,TIME,DT,NSTAGE,TORDER)
%
%      MASTER:                  Master structure
%      MESH:                    Mesh structure
%      UH(NPL,NC,NE):           Vector of unknowns (initial guess)
%      QH(NPL,NC,2,NE):         Vector of gradients of U (initial guess)
%      UHAT(NC,3*NPS,NF):       Vector of U_hat's (initial guess)
%      APP:                     Application structure
%      TIME:                    Time
%      DT:                      Timestep
%      NSTAGE:                  Number of stages
%      TORDER:                  Order of accuracy
%
%      UH(NPL,NC,NE):           Vector of unknowns
%      QH(NPL,NC,2,NE):         Vector of gradients of U
%      UHAT(NC,3*NPS,NF):       Vector of U_hat's
%
%      NPL:                     Number of DG nodes within an element
%      NC:                      Number of conservation equations solved (components)
%      NPS:                     Number of HDG nodes per edge (porder+1)
%      NE:                      Number of elements
%      NF:                      Number of faces
%
%  Ref: Diagonally Implicit Runge Kutta Methods for Stiff ODE's, 
%  by Roger Alexander, SINUM, vol. 14, no. 6, 1997.
%


ncu = pde.ncu;
[d,c,t] = dirkcoeff(nstage,torder);
d = d/dt;
Un = (1-sum(c))*UDG;
Pn = (1-sum(c))*PDG;
Hn = (1-sum(c))*UH;

% NOTE : for DIRK(3,3), c = [0,0,1] ; and for DIRK(2,2), c = [0,1].
% Therefore for these two schemes, (Un,Hn,Pn)=(0,0,0) at the beginning.

for istage = 1:nstage
    fprintf('DIRK stage :  %d\n', istage);
    fc_u = d(istage,istage);
    if pde.wave
        fc_q = fc_u;
    else
        fc_q = 1;
    end
    
    if istage==1
        SDG = d(1,1)*UDG; 
    elseif istage == 2
        SDG = (d(2,1)+d(2,2))*UDG - d(2,1)*UDGn;
    elseif istage == 3
        SDG = (d(3,1)+d(3,2)+d(3,3))*UDG - ...
              (d(3,1)/d(2,1))*((d(2,1)+d(2,2))*UDG-SDG) - d(3,2)*UDGn;    
    end        
    
    pde.time = time+dt*t(istage);
    pde.fc_u = fc_u;
    pde.fc_q = fc_q;
       
    [UDGn,UHn] = hdgsolve(master,mesh,pde,UDG,UH,SDG);
    Un = Un + c(istage)*UDGn;
    Hn = Hn + c(istage)*UHn;
    
    if pde.wave        
        if istage==1
            SPG = d(1,1)*PDG; 
        elseif istage == 2
            SPG = (d(2,1)+d(2,2))*PDG - d(2,1)*PDGn;
        elseif istage == 3
            SPG = (d(3,1)+d(3,2)+d(3,3))*PDG - ...
                  (d(3,1)/d(2,1))*((d(2,1)+d(2,2))*PDG-SPG) - d(3,2)*PDGn;    
        end                    
        PDGn = (1/d(istage,istage))*(UDGn(:,1:ncu,:)+SPG);
        Pn = Pn + c(istage)*PDGn;         
    end
end

function [d,c,t] = dirkcoeff(q,p)
% q - number of stages
% p - order of accuracy

% fprintf('DIRK (%d,%d) \n', p,q);

if q == 1 && p == 1
    a = 1;
    b = 1;
    t = 1;
elseif q == 1 && p == 2
    % fprintf('Warning: Innacutate for non-wave formuations - use (2,2) or (3,3)');
    a = 0.5;
    b = 1;
    t = 0.5;
elseif q == 2 && p == 2
    a = [1-0.5*sqrt(2), 0;
           0.5*sqrt(2), 1-0.5*sqrt(2)];
    b = [  0.5*sqrt(2), 1-0.5*sqrt(2)];
    t = [1-0.5*sqrt(2), 1];
elseif q == 2 && p == 3
    % fprintf('Warning: Innacutate for non-wave formuations - use (2,2) or (3,3)');
    a = [0.5+0.5/sqrt(3), 0;
              -1/sqrt(3), 0.5 + 0.5/sqrt(3)];
    b = [            0.5, 0.5];
    t = [0.5+0.5/sqrt(3), 0.5-0.5/sqrt(3)];
elseif q == 3 && p == 3
    a1 = 0.4358665215;
    t1 = (1+a1)/2;
    b1 = -(6*a1^2-16*a1+1)/4;
    b2 = (6*a1^2-20*a1+5)/4;    
    a = [a1   ,  0, 0;
         t1-a1, a1, 0;
         b1   , b2, a1];
    b = [b1, b2, a1];
    t = [a1, t1, 1];    
elseif q == 3 && p == 4
    % fprintf('Warning: Innacutate for non-wave formuations - use (2,2) or (3,3)');
    a1 = 2*cos(pi/18)/sqrt(3);
    a = [ 0.5*(1+a1),            0,          0;
             -0.5*a1,   0.5*(1+a1),          0;
                1+a1,    -(1+2*a1), 0.5*(1+a1)];
    b = [ 1/(6*a1^2), 1-1/(3*a1^2), 1/(6*a1^2)];
    t = [ 0.5*(1+a1),          0.5, 0.5*(1-a1)];
elseif q == 5 && p == 5
    a11 = (6-sqrt(6))/10;
    a21 = (-6+5*sqrt(6))/14;
    a31 = (888+607*sqrt(6))/2850; 
    a32 = (126-161*sqrt(6))/1425;
    a41 = (3153-3082*sqrt(6))/14250;
    a42 = (3213+1148*sqrt(6))/28500;
    a43 = (-267+88*sqrt(6))/500;
    a51 = (-32583+14638*sqrt(6))/71250;
    a52 = (-17199+364*sqrt(6))/142500;
    a53 = (1329-544*sqrt(6))/2500;
    a54 = (-96+131*sqrt(6))/625;
    b1  = 0;
    b2  = 0;
    b3  = (1/9);
    b4  = (16-sqrt(6))/36;
    b5  = (16+sqrt(6))/36;
    t1  = a11;
    t2  = (6+9*sqrt(6))/35;
    t3  = 1;
    t4  = (4-sqrt(6))/10;
    t5  = (4+sqrt(6))/10;
    a = [a11, 0, 0, 0, 0;
         a21, a11, 0, 0, 0;
         a31, a32, a11, 0, 0;
         a41, a42, a43, a11, 0;
         a51, a52, a53, a54, a11];
    b = [b1, b2, b3, b4, b5]; 
    t = [t1, t2, t3, t4, t5];            
else
    error('Invalid (q,p) combination');
end
d = inv(a);
c = b*d;

