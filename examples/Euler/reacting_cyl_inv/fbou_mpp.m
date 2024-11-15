function [fh,fh_udg,fh_uh] = fbou(ib,uinf,nl,p,udg,uh,param,time)
%FHAT flux function
%   [fh,fhu,fhq,fhm] = fhat(nl,p,u,q,m,param)
%
%      NL(N,ND)              Normal N points
%      P(N,ND)               Coordinates for N points
%      U(N,NC)               Unknown vector for N points with NC components
%      Q(N,NC,ND)            Flux vector for N points with NC components in the
%                            coordinate directions
%      M(N,NC)               Hybrid unkowns for N points with NC components
%      PARAM                 Parameter list
%      FH(N,NC):             Volume flux at N points
%      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
%      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
%      FHM(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Lambda

[ng,nch] = size(uh);
nc = size(udg,2);

% gam = param{1};
ns = 5;

u = udg(:,1:nch);

switch ib
    case 1  % Far field
        error("NOT IMPLEMENTED")
    case 2  % Wall
        un = u(:,ns+1).*nl(:,1) + u(:,ns+2).*nl(:,2);        
        uinf = u;
        uinf(:,ns+1) = uinf(:,ns+1) - nl(:,1).*un;
        uinf(:,ns+2) = uinf(:,ns+2) - nl(:,2).*un;
        
        fh = uinf - uh;
                
        fh_udg = zeros(ng,nch,nc);
        fh_udg(:,ns+1,ns+1:ns+2) = [ones(ng,1)-nl(:,1).*nl(:,1), -nl(:,1).*nl(:,2)]; 
        fh_udg(:,ns+2,ns+1:ns+2) = [ -nl(:,1).*nl(:,2), ones(ng,1)-nl(:,2).*nl(:,2)];
        fh_udg(:,ns+3,ns+3) = ones(ng,1); 
        for i = 1:5
            fh_udg(:,i,i) = ones(ng,1); 
        end
       
        fh_uh = zeros(ng,nch,nch);
        fh_uh(:,1,1) = -1;
        fh_uh(:,2,2) = -1;
        fh_uh(:,3,3) = -1;                
        fh_uh(:,4,4) = -1; 
        fh_uh(:,5,5) = -1; 
        fh_uh(:,6,6) = -1; 
        fh_uh(:,7,7) = -1; 
        fh_uh(:,8,8) = -1; 
    case 3 % Prescribed pressure
        error("NOT IMPLEMENTED")
    case 4 % Prescribed pressure
        error("NOT IMPLEMENTED")
    case 5                    
        fh = u - uh;
                
        fh_udg = zeros(ng,nch,nc);
        fh_udg(:,1,1) = ones(ng,1); 
        fh_udg(:,2,2) = ones(ng,1); 
        fh_udg(:,3,3) = ones(ng,1); 
        fh_udg(:,4,4) = ones(ng,1); 
        fh_udg(:,5,5) = ones(ng,1); 
        fh_udg(:,6,6) = ones(ng,1); 
        fh_udg(:,7,7) = ones(ng,1); 
        fh_udg(:,8,8) = ones(ng,1); 

       
        fh_uh = zeros(ng,nch,nch);
        fh_uh(:,1,1) = -1;
        fh_uh(:,2,2) = -1;
        fh_uh(:,3,3) = -1;                
        fh_uh(:,4,4) = -1;  
        fh_uh(:,5,5) = -1;
        fh_uh(:,6,6) = -1;
        fh_uh(:,7,7) = -1;                
        fh_uh(:,8,8) = -1;        
    case 6                    
        fh = uinf - uh;
                
        fh_udg = zeros(ng,nch,nc);
       
        fh_uh = zeros(ng,nch,nch);
        fh_uh(:,1,1) = -1;
        fh_uh(:,2,2) = -1;
        fh_uh(:,3,3) = -1;                
        fh_uh(:,4,4) = -1;     
        fh_uh(:,5,5) = -1;
        fh_uh(:,6,6) = -1;
        fh_uh(:,7,7) = -1;                
        fh_uh(:,8,8) = -1;               
    otherwise
        error('unknown boundary type');
end



