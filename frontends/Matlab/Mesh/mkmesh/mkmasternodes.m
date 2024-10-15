function [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface] = mkmasternodes(porder,dim,elemtype,nodetype)
%UNIFORMLOCALPNTS 2-D Mesh Generator for the master element.
%   [PLOCAL,TLOCAL]=UNIFORMLOCALPNTS(PORDER,NODETYPE)
%
%      PORDER    :  Orders of the complete polynomial 
%      DIM       :  Dimension
%      ELEMTYPE  :  Flag determining element type
%                   Flag = 0 tri/tet elements (default)
%                   Flag = 1 quad/hex elements
%      NODETYPE  :  Flag determining node distribution 
%                   Flag = 0 uniform distribution (default)
%                   Flag = 1 nonuniform distribution
%
%      PLOCAL    :  Node positions on the master volume element
%      TLOCAL    :  Element connectivity on the master vloume element
%      PLOCFC    :  Node positions on the master face element
%      PERMNODE  :  Node indices on the vertices of the master element
%      PERMEDFE  :  Node indices on the edges of the master element
%      PERMFACE  :  Node indices on the faces of the master element

if nargin<3, elemtype=0; end
if nargin<4, nodetype=0; end

switch dim
    case 1 % 1D        
        [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface] = nodes1d(porder(1),nodetype); 
    case 2 % 2D
        if elemtype==0     % tri        
            [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface] = trinodes2d(porder(1),nodetype);            
        elseif elemtype==1 % quad
            [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface] = quadnodes2d(porder,nodetype);         
        end        
    case 3 %3D
        if elemtype==0     % tet
            [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface] = tetnodes3d(porder(1),nodetype); 
        elseif elemtype==1 % hex
            [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface] = hexnodes3d(porder,nodetype);         
        end
    otherwise
        error('Only can handle dim=1, dim=2 or dim=3');
end


function [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface] = nodes1d(porder,nodetype)
% nodes on the interval [0, 1]

permedge = [];
permface = [];
plocfc   = [];
tlocfc   = [];
if porder==0
    plocal = 0.5;    
    tlocal = [1 2];
    permnode = 1;
    permedge = 1;
    permface = 1;
else    
    if nodetype==0
        plocal = linspace(0,1,porder+1)';    
    else
        plocal = xchenodes(porder);
    end    
    tlocal = [(1:porder)' (2:porder+1)'];  
    permnode = [1 porder+1];
    permface = permnode;
end


function [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface] = trinodes2d(porder,nodetype)
% nodes on the standard triangle 

permface = [];

if porder==0
    plocal = [1/3 1/3];
    tlocal = [1 2 3];
    plocfc = 1/2;
    tlocfc = [1 2];
    permnode = [1 1 1];      
    permedge = [1 1 1];
    permface = [1 1 1];
else    
    if nodetype==0
        [u,v]=ndgrid((0:porder)/porder,(0:porder)/porder);
        plocal=[u(:),v(:)];
        plocal=[1-sum(plocal,2),plocal];
        plocal=plocal(plocal(:,1)>=0,2:3);
    else 
        [plocal] = trinodes(porder);         
    end
    
    plocfc = plocal(1:porder+1,1);
    tlocfc = [(1:porder)' (2:porder+1)'];  
    
    tlocal=zeros(0,3);
    loc=0;
    for ii=1:porder
        jj=porder+1-ii;
        t1=[1:jj; 2:jj+1; jj+1+(1:jj)]';
        t2=[2:jj; jj+1+(2:jj); jj+1+(1:jj-1)]';
        tlocal=[tlocal; t1+loc; t2+loc];
        loc=loc+jj+1;
    end 
    
    permnode = [1 porder+1 0.5*(porder+1)*(porder+2)];      
    permedge(1,:) = [porder+1 0.5*(porder+1)*(porder+2) 1];
    for ii=2:porder+1
        permedge(ii,1) = permedge(ii-1,1) + porder+2-ii;  
        permedge(ii,2) = permedge(ii-1,2) - ii;
        permedge(ii,3) = ii;        
    end
end

function [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface] = tetnodes3d(porder,nodetype)
% nodes on the standard tetrahedron

if porder==0
    plocal=[1/4 1/4 1/4];
    tlocal=[1 2 3 4];
    plocfc=[1/3 1/3];
    tlocfc=[1 2 3];
    permnode = [1 1 1 1];      
    permedge = [1 1 1 1];    
    permface = [1 1 1 1];
else    
    if nodetype==0
        xi = nodes1d(porder,nodetype);    
        ploc2d = trinodes2d(porder,nodetype);
        plocal = [ploc2d zeros(size(ploc2d,1),1)];
        for ii=porder-1:-1:1
            ploc2d = trinodes2d(ii,nodetype);
            plocal = [plocal; [xi(ii+1)*ploc2d xi(porder-ii+1)*ones(size(ploc2d,1),1)]];
        end
        plocal = [plocal; [0 0 1]];    
    else
        plocal = tetnodes(porder);
    end
    
    plocfc = plocal(1:0.5*(porder+1)*(porder+2),1:2);
    tlocfc = zeros(0,3);
    loc=0;
    for ii=1:porder
        jj=porder+1-ii;
        t1=[1:jj; 2:jj+1; jj+1+(1:jj)]';
        t2=[2:jj; jj+1+(2:jj); jj+1+(1:jj-1)]';
        tlocfc=[tlocfc; t1+loc; t2+loc];
        loc=loc+jj+1;
    end
    
    [p,tlocal]=tetrahedronmesh(porder+1);
    
    permnode = [1 porder+1 0.5*(porder+1)*(porder+2) (porder+1)*(porder+2)*(porder+3)/6];      
    
    permedge(1,:) = [1 1 1 porder+1 0.5*(porder+1)*(porder+2) (porder+1)*(porder+2)*(porder+3)/6];
    for ii=2:porder+1
        permedge(ii,1) = ii;        
        permedge(ii,2) = permedge(ii-1,2) + porder+3-ii;  
        permedge(ii,3) = permedge(ii-1,3) + 0.5*(porder+3-ii)*(porder+4-ii);        
        permedge(ii,4) = permedge(ii-1,4) + porder+2-ii;  
        permedge(ii,5) = permedge(ii-1,5) + 0.5*(porder+2-ii)*(porder+3-ii);        
        permedge(ii,6) = permedge(ii-1,6) - 0.5*(ii-1)*(ii) - ii + 1;        
    end
        
    ploc = plocal;
    ind  = abs(ploc)<1e-7;
    ploc(ind) = 0;
    
    %face=[2,3,4]    
    ind = find(abs(ploc(:,1)+ploc(:,2)+ploc(:,3)-1)<1e-6);     
    permface(:,1)=ind;    
        
    %face=[1,4,3]
    ind = find(abs(ploc(:,1))<1e-6);         
    ii = faceindex(porder);
    permface(:,2)=ind(ii);        
    
    %face=[1,2,4]
    ind = find(abs(ploc(:,2))<1e-6);         
    permface(:,3)=ind;    
    
    %face=[1,3,2]
    ind = find(abs(ploc(:,3))<1e-6);               
    ii = faceindex(porder);
    permface(:,4)=ind(ii);    
end


function [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface] = prismnodes3d(porder,nodetype)
% nodes on the standard tetrahedron

if porder==0
    plocal=[1/3 1/3 1/2];
    tlocal=[1 2 3 4 5 6];
    plocfc=[1/3 1/3];
    tlocfc=[1 2 3];
    permnode = [1 1 1 1];      
    permedge = [1 1 1 1];    
    permface = [1 1 1 1];
else    
    if nodetype==0
        xi = nodes1d(porder,nodetype);    
        ploc2d = trinodes2d(porder,nodetype);
        plocal = [ploc2d zeros(size(ploc2d,1),1)];
        for ii=porder-1:-1:1
            ploc2d = trinodes2d(ii,nodetype);
            plocal = [plocal; [xi(ii+1)*ploc2d xi(porder-ii+1)*ones(size(ploc2d,1),1)]];
        end
        plocal = [plocal; [0 0 1]];    
    else
        plocal = tetnodes(porder);
    end
    
    plocfc = plocal(1:0.5*(porder+1)*(porder+2),1:2);
    tlocfc = zeros(0,3);
    loc=0;
    for ii=1:porder
        jj=porder+1-ii;
        t1=[1:jj; 2:jj+1; jj+1+(1:jj)]';
        t2=[2:jj; jj+1+(2:jj); jj+1+(1:jj-1)]';
        tlocfc=[tlocfc; t1+loc; t2+loc];
        loc=loc+jj+1;
    end
    
    [p,tlocal]=tetrahedronmesh(porder+1);
    
    permnode = [1 porder+1 0.5*(porder+1)*(porder+2) (porder+1)*(porder+2)*(porder+3)/6];      
    
    permedge(1,:) = [1 1 1 porder+1 0.5*(porder+1)*(porder+2) (porder+1)*(porder+2)*(porder+3)/6];
    for ii=2:porder+1
        permedge(ii,1) = ii;        
        permedge(ii,2) = permedge(ii-1,2) + porder+3-ii;  
        permedge(ii,3) = permedge(ii-1,3) + 0.5*(porder+3-ii)*(porder+4-ii);        
        permedge(ii,4) = permedge(ii-1,4) + porder+2-ii;  
        permedge(ii,5) = permedge(ii-1,5) + 0.5*(porder+2-ii)*(porder+3-ii);        
        permedge(ii,6) = permedge(ii-1,6) - 0.5*(ii-1)*(ii) - ii + 1;        
    end
        
    ploc = plocal;
    ind  = abs(ploc)<1e-7;
    ploc(ind) = 0;
    
    %face=[2,3,4]    
    ind = find(abs(ploc(:,1)+ploc(:,2)+ploc(:,3)-1)<1e-6);     
    permface(:,1)=ind;    
        
    %face=[1,4,3]
    ind = find(abs(ploc(:,1))<1e-6);         
    ii = faceindex(porder);
    permface(:,2)=ind(ii);        
    
    %face=[1,2,4]
    ind = find(abs(ploc(:,2))<1e-6);         
    permface(:,3)=ind;    
    
    %face=[1,3,2]
    ind = find(abs(ploc(:,3))<1e-6);               
    ii = faceindex(porder);
    permface(:,4)=ind(ii);    
end

function [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface] = quadnodes2d(porder,nodetype)
% nodes on the standard quad 

permface = [];

if porder==0
    plocal=[1/2 1/2];
    tlocal=[1 2 3 4];
    plocfc=[1/2];
    tlocfc=[1 2];
    permnode = [1 1 1 1];      
    permedge = [1 1 1 1];
    permface = [1 1 1 1];
else    
    if length(porder)==1
        porder = porder*ones(2,1);
    end
    xi  = nodes1d(porder(1),nodetype);
    eta  = nodes1d(porder(2),nodetype);
    [x,y]  = ndgrid(eta,xi);    
    plocal = [x(:) y(:)];
    
    plocfc = eta;
    tlocfc = [(1:porder(1))' (2:porder(1)+1)'];  
    
    m = porder(1)+1;
    n = porder(2)+1;
    t = [1 2 m+2 m+1];    
    t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
    t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');        
    tlocal = t;
    
    permnode = [1 m m*n m*(n-1)+1];          
    permedge(:,1) = 1:m;  
    permedge(:,2) = m:m:(m*n);
    permedge(:,3) = m*n:-1:(m*(n-1)+1);
    permedge(:,4) = (m*(n-1)+1):-m:1;
end


function [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface] = hexnodes3d(porder,nodetype)
% nodes on the standard hex

if porder==0
    plocal=[1/2 1/2 1/2];
    tlocal=[1 2 3 4 5 6 7 8];    
    plocfc=[1/2 1/2];
    tlocfc=[1 2 3 4];
    permnode = [1 1 1 1 1 1 1 1];      
    permedge = [1 1 1 1 1 1 1 1];   
    permface = [1 1 1 1 1 1 1 1];
else    
    if length(porder)==1
        porder = porder*ones(3,1);
    end
    xi  = nodes1d(porder(1),nodetype);
    eta  = nodes1d(porder(2),nodetype);
    zeta  = nodes1d(porder(3),nodetype);
    [x,y,z]  = ndgrid(xi,eta,zeta);    
    plocal = [x(:) y(:) z(:)];
    
    [x,y]  = ndgrid(xi,eta);    
    plocfc = [x(:) y(:)];
    
    m = porder(1)+1;
    n = porder(2)+1;
    t = [1 2 m+2 m+1];    
    t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
    t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');        
    tlocfc = t;
    
    m = porder(1)+1;
    n = porder(2)+1;
    k = porder(3)+1;
    t = [1 2 m+2 m+1 m*n+1 m*n+2 m*n+m+2 m*n+m+1];
    t=kron(t,ones(k-1,1))+kron(ones(size(t)),(0:k-2)'*(m*n));        
    t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
    t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');            
    tlocal = t;
    
    permnode = [1 m m*n m*(n-1)+1 m*n*(k-1)+1 m*n*(k-1)+m m*n*(k-1)+m*n m*n*(k-1)+m*(n-1)+1];          
    
    permedge(:,1) = 1:m;  
    permedge(:,2) = m:m:(m*n);
    permedge(:,3) = m*n:-1:(m*(n-1)+1);
    permedge(:,4) = (m*(n-1)+1):-m:1;
    permedge(:,5) = 1:m*n:(m*n*(k-1)+1);  
    permedge(:,6) = m:m*n:(m*n*(k-1)+m);
    permedge(:,7) = m*n:m*n:(m*n*k);
    permedge(:,8) = (m*(n-1)+1):m*n:(m*n*(k-1)+(m*(n-1)+1));
    permedge(:,9) = (m*n*(k-1)+1):(m*n*(k-1)+m);  
    permedge(:,10) = (m*n*(k-1)+m):m:(m*n*k);
    permedge(:,11) = m*n*k:-1:(m*n*k-m+1);
    permedge(:,12) = (m*n*k-m+1):-m:(m*n*(k-1)+1);
        
    % face=[1,4,3,2]
    f = krongrid((1:m:(m*(n-1)+1)), 0:1:m-1);
    permface(:,1) = f(:);    
    % face [5,6,7,8]
    f = krongrid((m*n*(k-1)+1):1:(m*n*(k-1)+m), 0:m:(m*(n-1)));
    permface(:,2) = f(:);
    % face [1,2,6,5]
    f = krongrid(1:m, 0:m*n:(m*n*(k-1)));
    permface(:,3) = f(:);
    % face [3,4,8,7]        
    f = krongrid(m*n:-1:(m*n-m+1), 0:m*n:(m*n*(k-1)));
    permface(:,4) = f(:);
    % face [2,3,7,6]
    f = krongrid(m:m:m*n, 0:m*n:(m*n*(k-1)));
    permface(:,5) = f(:);
    % face [4,1,5,8]
    f = krongrid((m*n-m+1):-m:1, 0:m*n:(m*n*(k-1)));
    permface(:,6) = f(:);
end

    
function xi = xchenodes(N)
% Extended Chebyshev nodes on the interval [0, 1]

if N==0
    xi = 0.5;
else
    n=N+1; k=1:n;       
    xi = -cos((2*k-1)*pi/(2*n))/cos(pi/(2*n));    
    xi = 0.5+0.5*xi';
end


function c = krongrid(a,b)

c = kron(a(:),ones(1,length(b))) + kron(ones(length(a),1),(b(:))');


function fi = faceindex(porder)

n = porder+1;
a = zeros(n,n);
k = 1;
for j=1:n
    for i=1:n+1-j        
        a(i,j) = k;
        k = k + 1;
    end
end
b = a';

in = b>0;
fi = b(in);
