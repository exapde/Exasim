function [p,t,p1]=cart2dg(elemtype,porder,X,Y,Z,opts)

if nargin==4
  opts=struct;
  dim=2;
elseif nargin==5 && isstruct(Z)
  opts=Z;
  dim=2;
elseif nargin==5
  opts=struct;
  dim=3;
elseif nargin==6
  dim=3;
end

if ~isfield(opts,'ind'), opts.ind=1; end

if elemtype==0
    if dim==2
      [m,n]=size(X);
      XX=permute(cat(3,X,Y),[3,1,2]);
      X1=XX(:,1:porder:m,1:porder:n);
      X1T=block2tets(X1,opts.ind);
      p1=reshape(X1T,2,[])';
      t1=reshape(1:numel(X1T)/2,3,size(X1T,3))';
      flip=simpvol(p1,t1)<0;
      t1(flip,[1,2])=t1(flip,[2,1]);

      [X0,Y0]=ndgrid(1:m,1:n);
      XX0=permute(cat(3,X0,Y0),[3,1,2]);
      X0=XX0(:,1:porder:m,1:porder:n);  
      X0T=block2tets(X0,opts.ind); 
      p0=reshape(X0T,2,[])';        
    elseif dim==3
      [m,n,o]=size(X);

      XX=permute(cat(4,X,Y,Z),[4,1,2,3]);
      X1=XX(:,1:porder:m,1:porder:n,1:porder:o);
      X1T=block2tets(X1,opts.ind);
      p1=reshape(X1T,3,[])';
      t1=reshape(1:numel(X1T)/3,4,size(X1T,3))';
      flip=simpvol(p1,t1)<0;
      t1(flip,[1,2])=t1(flip,[2,1]);

      [X0,Y0,Z0]=ndgrid(1:m,1:n,1:o);
      XX0=permute(cat(4,X0,Y0,Z0),[4,1,2,3]);
      X0=XX0(:,1:porder:m,1:porder:n,1:porder:o);
      X0T=block2tets(X0,opts.ind);
      p0=reshape(X0T,3,[])';
    end
elseif elemtype==1
    if dim==2
        [m,n]=size(X);
        x=X(1:porder:m,1:porder:n);
        y=Y(1:porder:m,1:porder:n);    
        [m,n]=size(x);
        p1 = [x(:) y(:)];
        %t1 = [1 2 m+2 m+1];
        t1 = [1 m+1 m+2 2];
        t1 = kron(t1,ones(n-1,1))+kron(ones(size(t1)),(0:n-2)'*m);
        t1 = kron(t1,ones(m-1,1))+kron(ones(size(t1)),(0:m-2)');
        
        [m,n]=size(X);        
        [X0,Y0]=ndgrid(1:m,1:n);        
        X0=X0(1:porder:m,1:porder:n);  
        Y0=Y0(1:porder:m,1:porder:n);          
        p0=[X0(:) Y0(:)];
    end        
    if dim==3
        [m,n,o]=size(X);
        x=X(1:porder:m,1:porder:n,1:porder:o);
        y=Y(1:porder:m,1:porder:n,1:porder:o);    
        z=Z(1:porder:m,1:porder:n,1:porder:o);    
        [m,n,o]=size(x);
        p1 = [x(:) y(:) z(:)];
        %t1 = [1 2 m+2 m+1 m*n+1 m*n+2 m*n+m+2 m*n+m+1];
        t1 = [1 m+1 m+2 2 m*n+1 m*n+m+1 m*n+m+2 m*n+2];
        t1=kron(t1,ones(o-1,1))+kron(ones(size(t1)),(0:o-2)'*(m*n));        
        t1=kron(t1,ones(n-1,1))+kron(ones(size(t1)),(0:n-2)'*m);
        t1=kron(t1,ones(m-1,1))+kron(ones(size(t1)),(0:m-2)');            
        
        [m,n,o]=size(X);        
        [X0,Y0,Z0]=ndgrid(1:m,1:n,1:o);        
        X0=X0(1:porder:m,1:porder:n,1:porder:o);  
        Y0=Y0(1:porder:m,1:porder:n,1:porder:o);
        Z0=Z0(1:porder:m,1:porder:n,1:porder:o);
        p0=[X0(:) Y0(:) Z0(:)];
    end        
end

s = localbasis(porder,dim,elemtype);
ns=size(s,1);

[nt1,nfe]=size(t1);
dgp0=zeros(ns,dim,nt1);
for comp=1:dim
  for trinode=1:nfe
    dp0=s(:,trinode)*p0(t1(:,trinode),comp)';
    dgp0(:,comp,:)=dgp0(:,comp,:)+permute(dp0,[1,3,2]);
  end
end

[m,n]=size(X);        
dgp0=round(dgp0);
if dim==2
  ix=sub2ind([m,n],dgp0(:,1,:),dgp0(:,2,:));
  dgp(:,1,:)=X(ix);
  dgp(:,2,:)=Y(ix);
elseif dim==3
  ix=sub2ind([m,n,o],dgp0(:,1,:),dgp0(:,2,:),dgp0(:,3,:));
  dgp(:,1,:)=X(ix);
  dgp(:,2,:)=Y(ix);
  dgp(:,3,:)=Z(ix);
end

if elemtype==0
    [p1,t1]=fixmesh(p1,t1);
else
    snap=max(max(p1,[],1)-min(p1,[],1),[],2)*1024*eps;
    [foo,ix,jx]=unique(round(p1/snap)*snap,'rows');
    p1=p1(ix,:);
    t1=jx(t1);
end
p=p1; t=t1; p1=dgp;

function philocal = localbasis(porder,nd,elemtype)

[plocal,tlocal] = masternodes(porder,nd,elemtype,0);
npv = size(tlocal,2);

if nd==1
    xi  = plocal(:,1);
    philocal(:,1) = 1 - xi;
    philocal(:,2) = xi;
elseif nd==2 && npv==3 % tri
    xi  = plocal(:,1);
    eta = plocal(:,2);    
    philocal(:,1) = 1 - xi - eta;
    philocal(:,2) = xi;
    philocal(:,3) = eta;
elseif nd==2 && npv==4 % quad
    xi  = plocal(:,1);
    eta = plocal(:,2);
    philocal(:,1) = (1-xi).*(1-eta);
    philocal(:,2) = xi.*(1-eta);
    philocal(:,3) = xi.*eta;
    philocal(:,4) = (1-xi).*eta;
elseif nd==3 && npv==4 % tet
    xi   = plocal(:,1);
    eta  = plocal(:,2);
    zeta = plocal(:,3);
    philocal(:,1) = 1 - xi - eta - zeta;
    philocal(:,2) = xi;
    philocal(:,3) = eta;
    philocal(:,4) = zeta;
elseif nd==3 && npv==8 % hex
    xi   = plocal(:,1);
    eta  = plocal(:,2);
    zeta = plocal(:,3);
    philocal(:,1) = (1-xi).*(1-eta).*(1-zeta);
    philocal(:,2) = xi.*(1-eta).*(1-zeta);
    philocal(:,3) = xi.*eta.*(1-zeta);
    philocal(:,4) = (1-xi).*eta.*(1-zeta);    
    philocal(:,5) = (1-xi).*(1-eta).*(zeta);
    philocal(:,6) = xi.*(1-eta).*(zeta);
    philocal(:,7) = xi.*eta.*(zeta);
    philocal(:,8) = (1-xi).*eta.*(zeta);        
end

