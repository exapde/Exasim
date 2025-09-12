function mastergrid = mkmastergriddual(mastergrid,pdual,ndual)


n = mastergrid.nref - ndual;

nd = mastergrid.nd;
porder   = mastergrid.porder;
elemtype = mastergrid.elemtype;
nodetype = mastergrid.nodetype;

mastergrid.ndual = ndual;
mastergrid.pdual = pdual;

% local nodes for pdual
mastergrid.plocvldual = masternodes(pdual,nd,elemtype,nodetype);

if n == 0
    % shape functions at the Gauss points of the basis set porder
    mastergrid.shapdual = mkshape(pdual,plocvldual,mastergrid.gpvl,mastergrid.elemtype);
elseif n>0
    % dual grid
    [mastergrid.pd,mastergrid.td] = masternodes(max(2^ndual,1),nd,elemtype,nodetype);
    
    % primal-dual grid
    [mastergrid.ps,mastergrid.ts] = masternodes(max(2^n,1),nd,elemtype,nodetype);
    
    % gauss nodes for the primal element
    mastergrid.gaussnodes = mknodes(mastergrid.ps,mastergrid.ts,mastergrid.gpvl);
    
    % shape functions of the standard master element at the gauss points
    for j=1:size(mastergrid.gaussnodes,3)
        tmp = mkshape(pdual,mastergrid.plocvldual,mastergrid.gaussnodes(:,:,j),elemtype);
        mastergrid.shapdual(:,:,j) = tmp(:,:,1);
    end
    
    % primal-dual map
    for i=1:size(mastergrid.t,1)
        t1 = mastergrid.t(i,:);
        p1 = mastergrid.p(t1,:);
        for j=1:size(mastergrid.td,1)
            td = mastergrid.td(j,:);
            pd = mastergrid.pd(td,:);
            for k=1:size(mastergrid.ts,1)
                ts = mastergrid.ts(k,:);
                ps = mastergrid.ps(ts,:);
                p2 = mapp(ps,pd);
                if norm(mean(p1)-mean(p2)) <= 1e-6
                    mastergrid.pdmap(i,1) = j;
                    mastergrid.pdmap(i,2) = k;                                            
                end
            end
        end
    end
else
    % dual grid
    [mastergrid.pd,mastergrid.td] = masternodes(max(2^ndual,1),nd,elemtype,nodetype);
    
    % primal-dual grid
    [mastergrid.ps,mastergrid.ts] = masternodes(max(2^abs(n),1),nd,elemtype,nodetype);
    
    % Gauss points and weights on the dual element  
    [mastergrid.gpvldual,mastergrid.gwvldual] = gaussquad(2*pdual,nd,elemtype);

    % gauss nodes for the dual element
    mastergrid.gaussnodes = mknodes(mastergrid.ps,mastergrid.ts,mastergrid.gpvldual);
    
    % shape functions of the standard master element at the gauss points
    for j=1:size(mastergrid.gaussnodes,3)        
        tmp = mkshape(porder,mastergrid.plocvl,mastergrid.gaussnodes(:,:,j),elemtype);
        mastergrid.shapdual(:,:,j) = tmp(:,:,1);
    end
   
    % primal-dual map
    for i=1:size(mastergrid.t,1)
        t1 = mastergrid.t(i,:);
        p1 = mastergrid.p(t1,:);
        for j=1:size(mastergrid.td,1)
            td = mastergrid.td(j,:);
            pd = mastergrid.pd(td,:);
            for k=1:size(mastergrid.ts,1)
                ts = mastergrid.ts(k,:);
                ps = mastergrid.ps(ts,:);
                p2 = mapp(ps,p1);            
                if norm(mean(pd)-mean(p2)) <= 1e-6
                    mastergrid.pdmap(i,k) = j;                    
                end
            end
        end
    end
end

function p=mapp(p,x)
%MAPP  ndim-linear map points
%

nx = size(x,1);
nd = size(x,2);

if nd==1
     r=[0,1];
     
     C=x(:,1)'/[1,1; r];
     
     p=[0*p+1,p]*C';
elseif (nd==2) && (nx==4) % 2-D bilinear map    
     r=[0,1,1,0];
     s=[0,0,1,1];
     
     C=x(:,1)'/[1,1,1,1; r; r.*s; s];
     D=x(:,2)'/[1,1,1,1; r; r.*s; s];

     px=p(:,1); py=p(:,2);
     p=[0*px+1, px, px.*py, py]*[C;D]';
elseif (nd==3) && (nx==8) % 3-D bilinear map    
     r=[0,1,1,0,0,1,1,0];
     s=[0,0,1,1,0,0,1,1];
     t=[0,0,0,0,1,1,1,1];
     
     C=x(:,1)'/[1,1,1,1,1,1,1,1; r; r.*s; s; t; r.*t; r.*s.*t; s.*t];
     D=x(:,2)'/[1,1,1,1,1,1,1,1; r; r.*s; s; t; r.*t; r.*s.*t; s.*t];
     E=x(:,3)'/[1,1,1,1,1,1,1,1; r; r.*s; s; t; r.*t; r.*s.*t; s.*t];
             
     px=p(:,1); py=p(:,2); pz=p(:,3);
     p=[0*px+1, px, px.*py, py, pz, px.*pz, px.*py.*pz, py.*pz]*[C;D;E]';
elseif (nd==2) && (nx==3) % 2-D affine map    
     r=[0,1,0];
     s=[0,0,1];
     
     C=x(:,1)'/[1,1,1; r; s];
     D=x(:,2)'/[1,1,1; r; s];

     px=p(:,1); py=p(:,2);
     p=[0*px+1, px, py]*[C;D]'; 
elseif (nd==3) && (nx==4) % 3-D affine map    
     r=[0,1,0,0];
     s=[0,0,1,0];
     t=[0,0,0,1];
     
     C=x(:,1)'/[1,1,1,1; r; s; t];
     D=x(:,2)'/[1,1,1,1; r; s; t];
     E=x(:,3)'/[1,1,1,1; r; s; t];

     px=p(:,1); py=p(:,2); pz=p(:,3);
     p=[0*px+1, px, py, pz]*[C;D;E]';
end






