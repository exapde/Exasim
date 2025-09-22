function [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface,face] = mkmasternodes(porder,dim,elemtype,nodetype)
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
        [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface,face] = nodes1d(porder(1),nodetype); 
    case 2 % 2D
        if elemtype==0     % tri        
            [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface,face] = trinodes2d(porder(1),nodetype);            
        elseif elemtype==1 % quad
            [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface,face] = quadnodes2d(porder,nodetype);         
        end        
    case 3 %3D
        if elemtype==0     % tet
            [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface,face] = tetnodes3d(porder(1),nodetype); 
        elseif elemtype==1 % hex
            [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface,face] = hexnodes3d(porder,nodetype);         
        elseif elemtype==2 % prism
            [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface,face] = prismnodes3d(porder,nodetype);                     
        elseif elemtype==3 % pyramid
            [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface,face] = pyramidnodes3d(porder,nodetype);                                 
        end
    otherwise
        error('Only can handle dim=1, dim=2 or dim=3');
end


function [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface,face] = nodes1d(porder,nodetype)
% nodes on the interval [0, 1]

face{1}(1) = 1;
plocfc{1}(1) = 0;
tlocfc{1}(1) = 1;
face{2}(1) = 2;
plocfc{2}(1) = 0;
tlocfc{2}(1) = 1;
if porder==0
    plocal = 0.5;    
    tlocal = [1 2];
    permnode = 1;    
    permface{1}(1) = 1;
    permface{2}(1) = 1;
else    
    if nodetype==0
        plocal = linspace(0,1,porder+1)';    
    else
        plocal = xchenodes(porder);
    end    
    tlocal = [(1:porder)' (2:porder+1)'];  
    permnode = [1 porder+1]';
    permface{1}(1) = 1;
    permface{2}(1) = porder+1;
end
permedge = permnode;

function [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface,face] = trinodes2d(porder,nodetype)
% nodes on the standard triangle 

face{1} = [2 3];
face{2} = [3 1];
face{3} = [1 2];
if porder==0
    plocal = [1/3 1/3];
    tlocal = [1 2 3];
    plocfc = 1/2;
    tlocfc = [1 2];
    permnode = [1 1 1];      
    permedge = [1 1 1];      
    permface{1}(1) = 1;
    permface{2}(1) = 1;
    permface{3}(1) = 1;
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
    permface{1} = zeros(porder+1,1);
    permface{2} = zeros(porder+1,1);
    permface{3} = zeros(porder+1,1);
    permface{1}(1) = porder+1;
    permface{2}(1) = 0.5*(porder+1)*(porder+2);
    permface{3}(1) = 1;
    for ii=2:porder+1
        permface{1}(ii) = permface{1}(ii-1) + porder+2-ii;  
        permface{2}(ii) = permface{2}(ii-1) - ii;
        permface{3}(ii) = ii;        
    end
    permedge = [permface{1} permface{2} permface{3}];
end
for i = 1:3
    pface{i} = plocfc;
    tface{i} = tlocfc;
end
plocfc = pface;
tlocfc = tface;


function [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface,face] = tetnodes3d(porder,nodetype)
% nodes on the standard tetrahedron

% [[2,3,4];[1,4,3];[1,2,4];[1,3,2]] 
face{1} = [2 3 4];
face{2} = [1 4 3];
face{3} = [1 2 4];
face{4} = [1 3 2];
if porder==0
    plocal=[1/4 1/4 1/4];
    tlocal=[1 2 3 4];
    plocfc=[1/3 1/3];
    tlocfc=[1 2 3];
    permnode = [1 1 1 1];      
    permedge = [1 1 1 1];    
    permface{1} = 1;
    permface{2} = 1;
    permface{3} = 1;
    permface{4} = 1;
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
    
    [~,tlocal]=tetrahedronmesh(porder+1);
    
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
    permface{1}=ind;    
        
    %face=[1,4,3]
    ind = find(abs(ploc(:,1))<1e-6);         
    ii = faceindex(porder);
    permface{2}=ind(ii);        
    
    %face=[1,2,4]
    ind = find(abs(ploc(:,2))<1e-6);         
    permface{3}=ind;    
    
    %face=[1,3,2]
    ind = find(abs(ploc(:,3))<1e-6);               
    ii = faceindex(porder);
    permface{4}=ind(ii);    
end
for i = 1:4
    pface{i} = plocfc;
    tface{i} = tlocfc;
end
plocfc = pface;
tlocfc = tface;


function [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface,face] = prismnodes3d(porder,nodetype)
% nodes on the standard tetrahedron

%face=[[0,2,1];[3,4,5];[1,2,5,4];[2,0,3,5];[0,1,4,3]]+1
face{1} = [0,2,1]+1;
face{2} = [3,4,5]+1;
face{3} = [1,2,5,4]+1;
face{4} = [2,0,3,5]+1;
face{5} = [0,1,4,3]+1;
if sum(porder)==0
    plocal=[1/3 1/3 1/2];
    tlocal=[1 2 3 4 5 6];    
    plocfc{1}=[1/3 1/3];
    tlocfc{1}=[1 2 3];
    plocfc{2}=[1/3 1/3];
    tlocfc{2}=[1 2 3];
    plocfc{3}=[1/2 1/2];
    tlocfc{3}=[1 2 3 4];
    plocfc{4} = plocfc{3};    
    tlocfc{4} = tlocfc{3};
    plocfc{5} = plocfc{3};    
    tlocfc{5} = tlocfc{3};    
    permnode = [1 1 1 1 1 1];      
    permedge = [1 1 1 1 1 1 1 1 1];    
    for i=1:5                
        permface{i} = 1;
    end
else    
    if length(porder)==1
        porder = porder*ones(2,1);
    end
    [ploc2d,tloc2d,~,~,~,permedgetri] = trinodes2d(porder(1),nodetype);
    xi  = nodes1d(porder(1),nodetype); 
    zeta  = nodes1d(porder(2),nodetype);    
    plocal = [ploc2d zeta(1)*ones(size(ploc2d,1),1)];    
    for i = 2:length(zeta)
        plocal = [plocal; [ploc2d zeta(i)*ones(size(ploc2d,1),1)]];
    end    
    tlocal = zeros(size(tloc2d,1)*porder(2),6);
    
    plocfc{1} = plocal(1:0.5*(porder(1)+1)*(porder(1)+2),1:2);
    tlocfc{1} = zeros(0,3);
    loc=0;
    for ii=1:porder(1)
        jj=porder(1)+1-ii;
        t1=[1:jj; 2:jj+1; jj+1+(1:jj)]';
        t2=[2:jj; jj+1+(2:jj); jj+1+(1:jj-1)]';
        tlocfc{1}=[tlocfc{1}; t1+loc; t2+loc];
        loc=loc+jj+1;
    end
    plocfc{2} = plocfc{1};
    tlocfc{2} = tlocfc{1};
    
    [x,y]  = ndgrid(xi,zeta);    
    plocfc{3} = [x(:) y(:)];    
    m = length(xi);
    k = length(zeta);
    t = [1 2 m+2 m+1];    
    t=kron(t,ones(k-1,1))+kron(ones(size(t)),(0:k-2)'*m);
    t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');        
    tlocfc{3} = t;            
    plocfc{4} = plocfc{3};    
    tlocfc{4} = tlocfc{3};
    plocfc{5} = plocfc{3};    
    tlocfc{5} = tlocfc{3};
    
    n1z = porder(2)+1;
    n1d = porder(1)+1;
    n2d = 0.5*(porder(1)+1)*(porder(1)+2);    
    permnode = [1 n1d n2d n2d*(n1z-1)+1 n2d*(n1z-1)+n1d n2d*(n1z-1)+n2d];        
    
    permedge = zeros(max([n1d n1z]),9);
    permedge(1:n1d,1:3) = permedgetri;
    permedge(1:n1d,7:9) = permedgetri+n2d*(n1z-1);    
    for i=1:n1z            
        permedge(i,4:6) = [1 n1d n2d]+(i-1)*n2d;
    end
    
    ploc = plocal;
    ind  = abs(ploc)<1e-7;
    ploc(ind) = 0;
    
    %face=[1,3,2]
    ind = find(abs(ploc(:,3))<1e-6);               
    ii = faceindex(porder(1));
    permface{1}=ind(ii);    
    
    %face=[4,5,6]
    ind = find(abs(ploc(:,3))>1-1e-6);                   
    permface{2}=ind;    
    
    %face=[2,3,6,5]
    ind = find(abs(ploc(:,1)+ploc(:,2))>1-1e-6);                   
    permface{3}=ind;    
    
    %face=[3,1,4,6]
    f = krongrid(permedge(1:n1d,2), 0:n2d:(n2d*(n1z-1)));
    permface{4} = f(:);    
    
    %face=[1,2,5,4]
    ind = find(abs(ploc(:,2))<1e-6);                   
    permface{5}=ind;        
end

function [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface,face] = pyramidnodes3d(porder,nodetype)
% nodes on the standard tetrahedron

%face=[[0,3,2,1];[0,1,4];[1,2,4];[2,3,4];[3,0,4]]+1                
face{1} = [1,4,3,2];
face{2} = [1,2,5];
face{3} = [2,3,5];
face{4} = [3,4,5];
face{5} = [4,1,5];
if sum(porder)==0
    plocal=[1/2 1/2 1/2];
    tlocal=[1 2 3 4 5];
    plocfc{1}=[1/2 1/2];
    tlocfc{1}=[1 2 3 4];
    plocfc{2}=[1/3 1/3];
    tlocfc{2}=[1 2 3];
    plocfc{3} = plocfc{2};    
    tlocfc{3} = tlocfc{2};
    plocfc{4} = plocfc{2};    
    tlocfc{4} = tlocfc{2};
    plocfc{5} = plocfc{2};    
    tlocfc{5} = tlocfc{2};    
    permnode = [1 1 1 1 1];      
    permedge = [1 1 1 1 1 1 1 1];    
    for i=1:5        
        permface{i} = 1;
    end
else
    numNodes1d = porder + 1;
    nfe = 5;
    npv = (porder+1)*(porder+2)*(2*porder+3)/6;    
    
    [ptri,ttri] = trinodes2d(porder,nodetype);
    [pquad,tquad] = quadnodes2d(porder,nodetype);
    plocfc{1} = pquad;
    tlocfc{1} = tquad;
    for i=2:nfe
        plocfc{i} = ptri;
        tlocfc{i} = ttri;
    end
        
    plocal = zeros(npv,3);
    p1dz = nodes1d(porder,nodetype);
    m = 0;
    for i=1:numNodes1d
        if i<numNodes1d
            a = (i-1)/porder;
            x1 = a*0.5;
            x2 = 1 - a*0.5;
            p1dx = nodes1d(numNodes1d-i,nodetype);
            p1dx = x1 + (x2-x1)*p1dx;
        else
            p1dx = 0.5;
        end                
        for j = 1:(numNodes1d+1-i)
            for k = 1:(numNodes1d+1-i)
                plocal(m+(j-1)*(numNodes1d+1-i)+k,1) = p1dx(k);
                plocal(m+(j-1)*(numNodes1d+1-i)+k,2) = p1dx(j);
                plocal(m+(j-1)*(numNodes1d+1-i)+k,3) = p1dz(i);
            end
        end
        m = m + (numNodes1d+1-i)*(numNodes1d+1-i);                
    end
    tlocal = [1 2 3 4 5];          
            
    n1d = porder(1)+1;
    n2d = 0.5*(porder(1)+1)*(porder(1)+2);    
    n2q = (porder + 1) * (porder + 1);
    permnode = [1 n1d n2q n2q-porder npv];        
        
    permedge = zeros(n1d,8);
    m = n1d; n = n1d;
    permedge(1:m,1) = 1:m;  
    permedge(1:n,2) = m:m:(m*n);
    permedge(1:m,3) = m*n:-1:(m*(n-1)+1);
    permedge(1:n,4) = (m*(n-1)+1):-m:1;    
    m = 0;
    for i=1:n1d                    
        n = (n1d+1-i);        
        permedge(i,5:8) = [1 n n*n n*n-n+1]+m;
        m = m + n*n;        
    end
    
    % face=[1,4,3,2]
    m = n1d; n = n1d;
    f = krongrid((1:m:(m*(n-1)+1)), 0:1:m-1);
    permface{1} = f(:);    
    
    ploc = plocal;
    ind  = abs(ploc)<1e-7;
    ploc(ind) = 0;
    
    %face=[1,2,5]
    ind = find(abs(ploc(:,3)-2*ploc(:,2))<1e-6);       
    permface{2}=ind;    
    
    %face=[2,3,5]
    ind = find(abs(ploc(:,1)+0.5*ploc(:,3)-1)<1e-6);        
    permface{3}=ind;    
    
    %face=[3,4,5]
    ind = find(abs(ploc(:,2)+0.5*ploc(:,3)-1)<1e-6);          
    ii = faceindex(porder);        
    permface{4}=ind(ii);    
    
    %face=[4,1,5]
    ind = find(abs(ploc(:,3)-2*ploc(:,1))<1e-6);                   
    permface{5}=ind(ii);        
end


function [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface,face] = quadnodes2d(porder,nodetype)
% nodes on the standard quad 

% [[1,2];[2,3];[3,4];[4,1]]
face{1} = [1 2];
face{2} = [2 3];
face{3} = [3 4];
face{4} = [4 1];
if sum(porder)==0
    plocal=[1/2 1/2];
    tlocal=[1 2 3 4];    
    permnode = [1 1 1 1];      
    permedge = [1 1 1 1];
    for i = 1:4
        plocfc{i}=1/2;
        tlocfc{i}=[1 2];
        permface{i} = 1;
    end    
else    
    if length(porder)==1
        porder = porder*ones(2,1);
    end
    xi  = nodes1d(porder(1),nodetype);
    eta  = nodes1d(porder(2),nodetype);
    [x,y]  = ndgrid(xi,eta);    
    plocal = [x(:) y(:)];
        
    m = porder(1)+1;
    n = porder(2)+1;
    t = [1 2 m+2 m+1];    
    t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
    t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');        
    tlocal = t;
    
    plocfc{1} = xi;
    tlocfc{1} = [(1:porder(1))' (2:porder(1)+1)'];      
    plocfc{2} = eta;
    tlocfc{2} = [(1:porder(2))' (2:porder(2)+1)'];  
    plocfc{3} = xi;
    tlocfc{3} = [(1:porder(1))' (2:porder(1)+1)'];      
    plocfc{4} = eta;
    tlocfc{4} = [(1:porder(2))' (2:porder(2)+1)'];  
    
    permnode = [1 m m*n m*(n-1)+1];          
    permedge = zeros(max([m,n]),4);
    permedge(1:m,1) = 1:m;  
    permedge(1:n,2) = m:m:(m*n);
    permedge(1:m,3) = m*n:-1:(m*(n-1)+1);
    permedge(1:n,4) = (m*(n-1)+1):-m:1;
    permface{1} = permedge(1:m,1);
    permface{2} = permedge(1:n,2);
    permface{3} = permedge(1:m,3);
    permface{4} = permedge(1:n,4);
end


function [plocal,tlocal,plocfc,tlocfc,permnode,permedge,permface,face] = hexnodes3d(porder,nodetype)
% nodes on the standard hex

%face=[[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]]
face{1} = [1,4,3,2];
face{2} = [5,6,7,8];
face{3} = [1,2,6,5];
face{4} = [3,4,8,7];
face{5} = [2,3,7,6];
face{6} = [4,1,5,8];
if sum(porder)==0
    plocal=[1/2 1/2 1/2];
    tlocal=[1 2 3 4 5 6 7 8];        
    permnode = [1 1 1 1 1 1 1 1];      
    permedge = [1 1 1 1 1 1 1 1 1 1 1 1];   
    for i = 1:6
        plocfc{i}=[1/2 1/2];
        tlocfc{i}=[1 2 3 4];
        permface{i} = 1;
    end    
else    
    if length(porder)==1
        porder = porder*ones(3,1);
    end
    xi  = nodes1d(porder(1),nodetype);
    eta  = nodes1d(porder(2),nodetype);
    zeta  = nodes1d(porder(3),nodetype);
    [x,y,z]  = ndgrid(xi,eta,zeta);    
    plocal = [x(:) y(:) z(:)];
        
    m = porder(1)+1;
    n = porder(2)+1;
    k = porder(3)+1;
    t = [1 2 m+2 m+1 m*n+1 m*n+2 m*n+m+2 m*n+m+1];
    t=kron(t,ones(k-1,1))+kron(ones(size(t)),(0:k-2)'*(m*n));        
    t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
    t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');            
    tlocal = t;
        
    [x,y]  = ndgrid(eta,xi);    
    plocfc{1} = [x(:) y(:)];    
    t = [1 2 n+2 n+1];    
    t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)'*n);
    t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)');        
    tlocfc{1} = t;    
    
    [x,y]  = ndgrid(xi,eta);    
    plocfc{2} = [x(:) y(:)];    
    t = [1 2 m+2 m+1];    
    t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
    t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');        
    tlocfc{2} = t;    
    
    [x,y]  = ndgrid(xi,zeta);    
    plocfc{3} = [x(:) y(:)];    
    t = [1 2 m+2 m+1];    
    t=kron(t,ones(k-1,1))+kron(ones(size(t)),(0:k-2)'*m);
    t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');        
    tlocfc{3} = t;            
    plocfc{4} = plocfc{3};    
    tlocfc{4} = tlocfc{3};
        
    [x,y]  = ndgrid(eta,zeta);    
    plocfc{5} = [x(:) y(:)];    
    t = [1 2 n+2 n+1];    
    t=kron(t,ones(k-1,1))+kron(ones(size(t)),(0:k-2)'*n);
    t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)');        
    tlocfc{5} = t;    
    plocfc{6} = plocfc{5};    
    tlocfc{6} = tlocfc{5};
    
    permnode = [1 m m*n m*(n-1)+1 m*n*(k-1)+1 m*n*(k-1)+m m*n*(k-1)+m*n m*n*(k-1)+m*(n-1)+1];          
    
    permedge = zeros(max([m,n,k]),12);
    permedge(1:m,1) = 1:m;  
    permedge(1:n,2) = m:m:(m*n);
    permedge(1:m,3) = m*n:-1:(m*(n-1)+1);
    permedge(1:n,4) = (m*(n-1)+1):-m:1;
    permedge(1:k,5) = 1:m*n:(m*n*(k-1)+1);  
    permedge(1:k,6) = m:m*n:(m*n*(k-1)+m);
    permedge(1:k,7) = m*n:m*n:(m*n*k);
    permedge(1:k,8) = (m*(n-1)+1):m*n:(m*n*(k-1)+(m*(n-1)+1));
    permedge(1:m,9) = (m*n*(k-1)+1):(m*n*(k-1)+m);  
    permedge(1:n,10) = (m*n*(k-1)+m):m:(m*n*k);
    permedge(1:m,11) = m*n*k:-1:(m*n*k-m+1);
    permedge(1:n,12) = (m*n*k-m+1):-m:(m*n*(k-1)+1);
        
    % face=[1,4,3,2]
    f = krongrid((1:m:(m*(n-1)+1)), 0:1:m-1);
    permface{1} = f(:);    
    % face [5,6,7,8]
    f = krongrid((m*n*(k-1)+1):1:(m*n*(k-1)+m), 0:m:(m*(n-1)));
    permface{2} = f(:);
    % face [1,2,6,5]
    f = krongrid(1:m, 0:m*n:(m*n*(k-1)));
    permface{3} = f(:);
    % face [3,4,8,7]        
    f = krongrid(m*n:-1:(m*n-m+1), 0:m*n:(m*n*(k-1)));
    permface{4} = f(:);
    % face [2,3,7,6]
    f = krongrid(m:m:m*n, 0:m*n:(m*n*(k-1)));
    permface{5} = f(:);
    % face [4,1,5,8]
    f = krongrid((m*n-m+1):-m:1, 0:m*n:(m*n*(k-1)));
    permface{6} = f(:);
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







