function [f,fx,fy,fz] = koornwinder(x,porder)
%KOORNWINDER Vandermonde matrices for Koornwinder polynomials on simplex elements. 
%            In 1D the Legendre polynomilas normalized to [0,1] are used. 
%            In 2D the koornwinder basis is normalized to the (0,0),(1,0),(0,1) triangle. 
%            In 3D the koornwinder basis is normalized to the 
%            (0,0,0),(1,0,0),(0,1,0),(0,0,1) tetrahedron. The basis are orthonormal.
%
%   [F,FX,FY,FZ]=KOORNWINDER(X,P)
%
%      X:         Coordinates of the points wherethe polynomials 
%                 are to be evaluated (npoints,dim)
%      PORDER:    Maximum order of the polynomials consider. That
%                 is all polynomials of complete degree up to p,
%                 npoly = (PORDER+1)*(PORDER+2)/2
%      F:         Vandermonde matrix (npoints,npoly)
%      FX:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. x (npoints,npoly)
%      FY:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. y (npoints,npoly)
%      FZ:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. z (npoints,npoly)
%

dim=size(x,2);

switch dim
 case 1
     [f,fx]=koornwinder1d(x,porder);      % 1D 
     fy=[];
     fz=[];
 case 2
     [f,fx,fy]=koornwinder2d(x,porder);   % 2D
     fz=[];
 case 3
     [f,fx,fy,fz]=koornwinder3d(x,porder);% 3D           
 otherwise
     error('Only can handle dim=1, dim=2 or dim=3');
end


function [f,fx] = koornwinder1d(x,p)
%KOORNWINDER1D Vandermonde matrix for Legenedre polynomials in [0,1]
%   [F,FX]=KOORNWINDER(X,P)
%
%      X:         Coordinates of the points wherethe polynomials 
%                 are to be evaluated (NPOINTS)
%      P:         Maximum order of the polynomials consider. That
%                 is all polynomials of degree up to P, NPOLY=P+1
%      F:         Vandermonde matrix (NPOINTS,NPOLY)
%      FX:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. x (NPOINTS,NPOLY)
%

x = 2*x-1;

%pp=p+1;
f=zeros(numel(x),p+1);
fx=zeros(numel(x),p+1);

for ii=0:p
    pp = jacobi(ii,0,0);
% Normalization factor to ensure integration to one
    pp = pp*sqrt(2*ii+1);  
    dpp = polyder(pp); 
    pval  = polyval(pp,x);
    dpval = polyval(dpp,x);
    f(:,ii+1) = pval;
    fx(:,ii+1) = dpval;
end

fx = 2*fx;


function [f,fx,fy] = koornwinder2d(x,p)
%KOORNWINDER2D Vandermonde matrix for Koornwinder polynomials in 
%              the master triangle [0,0]-[1,0]-[0,1]
%   [F,FX,FY]=KOORNWINDER(X,P)
%
%      X:         Coordinates of the points wherethe polynomials 
%                 are to be evaluated (npoints,dim)
%      PORDER:    Maximum order of the polynomials consider. That
%                 is all polynomials of complete degree up to p,
%                 npoly = (PORDER+1)*(PORDER+2)/2
%      F:         Vandermonde matrix (npoints,npoly)
%      FX:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. x (npoints,npoly)
%      FY:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. y (npoints,npoly)
%
%
x = 2*x-1;
npol=prod((double(p)+(1:2))./(1:2));
f=zeros(size(x,1),npol);
fx=zeros(size(x,1),npol);
fy=zeros(size(x,1),npol);

pq = pascalindex2d(npol);

xc = x;
xc(:,2) = min( 0.99999999, xc(:,2)); % avoid singularity

e(:,1) = 2*(1+xc(:,1))./(1-xc(:,2))-1;
e(:,2) = xc(:,2);
ii = find(x(:,2) == 1);
% Correct values for function evaluation
e(ii,1) = -1;
e(ii,2) =  1;

for ii=1:npol
    pp = jacobi(pq(ii,1),0,0);
    qp = jacobi(pq(ii,2),2*pq(ii,1)+1,0);
    for i=1:pq(ii,1)
        qp = conv([-0.5,0.5],qp);
    end
    
    pval = polyval(pp,e(:,1));
    qval = polyval(qp,e(:,2));
    
    % Normalization factor to ensure integration to one    
    fc = sqrt((2*pq(ii,1)+1)*2*(pq(ii,1)+pq(ii,2)+1));

    f(:,ii) = fc*pval.*qval;
end


% Use displaced coordinate for derivative evaluation
e(:,1) = 2*(1+xc(:,1))./(1-xc(:,2))-1;
e(:,2) = xc(:,2);
de1(:,1) = 2./(1-xc(:,2));
de1(:,2) = 2*(1+xc(:,1))./(1-xc(:,2)).^2;

for ii=1:npol
    pp = jacobi(pq(ii,1),0,0);
    qp = jacobi(pq(ii,2),2*pq(ii,1)+1,0);
    for i=1:pq(ii,1)
        qp = conv([-0.5,0.5],qp);
    end

    dpp = polyder(pp);
    dqp = polyder(qp);

    pval = polyval(pp,e(:,1));
    qval = polyval(qp,e(:,2));
    
    dpval = polyval(dpp,e(:,1));
    dqval = polyval(dqp,e(:,2));
    
    % Normalization factor to ensure integration to one    
    fc = sqrt((2*pq(ii,1)+1)*2*(pq(ii,1)+pq(ii,2)+1));

    fx(:,ii) = fc*dpval.*qval.*de1(:,1);
    fy(:,ii) = fc*(dpval.*qval.*de1(:,2) + pval.*dqval);
end
fx = 2*fx;
fy = 2*fy;

function pq = pascalindex2d(p)
l=1;
for i=0:p
    for j=0:i
        pq(l,1)=i-j;
        pq(l,2)=j;
        l = l+1;
        if l>p 
           return;
        end
    end
end

function [f,fx,fy,fz] = koornwinder3d(x,p)
%KOORNWINDER2D Vandermonde matrix for Koornwinder polynomials in 
%              the master tetrahedron [0,0,0]-[1,0,0]-[0,1,0]-[0,0,1]
%   [F,FX,FY]=KOORNWINDER(X,P)
%
%      X:         Coordinates of the points where the polynomials 
%                 are to be evaluated (npoints,dim)
%      PORDER:    Maximum order of the polynomials consider. That
%                 is all polynomials of complete degree up to p,
%                 npoly = (PORDER+1)*(PORDER+2)*(PORDER+3)/6
%      F:         Vandermonde matrix (npoints,npoly)
%      FX:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. x (npoints,npoly)
%      FY:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. y (npoints,npoly)
%      FZ:        Vandermonde matrix for the derivative of the Koornwinder
%                 polynomials w.r.t. Z (npoints,npoly)

%
x = 2*x-1;
npol=(p+1)*(p+2)*(p+3)/6;
f=zeros(size(x,1),npol);
fx=zeros(size(x,1),npol);
fy=zeros(size(x,1),npol);
fz=zeros(size(x,1),npol);

pq = pascalindex3d(npol);

e  = x;
for ii=1:size(x,1)
    if (x(ii,2)+x(ii,3) ~= 0)
        e(ii,1) = -2*(1+x(ii,1))./(x(ii,2)+x(ii,3))-1;
    else
        e(ii,1) = -1;
    end
    if x(ii,3) ~= 1
        e(ii,2) = 2*(1+x(ii,2))./(1-x(ii,3))-1;
    else
        e(ii,2) = -1;
    end        
    e(ii,3) = x(ii,3);    
end

for ii=1:npol
    pp = jacobi(pq(ii,1),0,0);
    qp = jacobi(pq(ii,2),2*pq(ii,1)+1,0);
    rp = jacobi(pq(ii,3),2*pq(ii,1)+2*pq(ii,2)+2,0);
    for i=1:pq(ii,1)
        qp = conv([-0.5,0.5],qp);
    end
    for i=1:pq(ii,1)+pq(ii,2)        
        rp = conv([-0.5,0.5],rp);        
    end
    
    pval = polyval(pp,e(:,1));
    qval = polyval(qp,e(:,2));
    rval = polyval(rp,e(:,3));
    
    % Normalization factor to ensure integration to one    
    fc = sqrt((2*pq(ii,1)+1)*2*(pq(ii,1)+pq(ii,2)+1)*2*(pq(ii,1)+pq(ii,2)+pq(ii,3)+2));
    
    f(:,ii) = fc*pval.*qval.*rval;
end

% Use displaced coordinate for derivative evaluation
xc = x;
for ii=1:size(x,1)
    if (x(ii,2)+x(ii,3) == 0)        
        xc(ii,3) = -1e-8-xc(ii,2);
    end
    if x(ii,3) == 1
        xc(ii,3) = 0.99999999;
    end                
end

e(:,1) = -2*(1+xc(:,1))./(xc(:,2)+xc(:,3))-1;
e(:,2) = 2*(1+xc(:,2))./(1-xc(:,3))-1;
e(:,3) = xc(:,3);
de1(:,1) = -2./(xc(:,2)+xc(:,3));
de1(:,2) = 2*(1+xc(:,1))./(xc(:,2)+xc(:,3)).^2;
de1(:,3) = 2*(1+xc(:,1))./(xc(:,2)+xc(:,3)).^2;
de2(:,1) = 0*xc(:,1); 
de2(:,2) = 2./(1-xc(:,3));
de2(:,3) = 2*(1+xc(:,2))./(1-xc(:,3)).^2;

for ii=1:npol
    pp = jacobi(pq(ii,1),0,0);
    qp = jacobi(pq(ii,2),2*pq(ii,1)+1,0);
    rp = jacobi(pq(ii,3),2*pq(ii,1)+2*pq(ii,2)+2,0);
    for i=1:pq(ii,1)
        qp = conv([-0.5,0.5],qp);
    end
    for i=1:pq(ii,1)+pq(ii,2)        
        rp = conv([-0.5,0.5],rp);        
    end

    dpp = polyder(pp);
    dqp = polyder(qp);
    drp = polyder(rp);

    pval = polyval(pp,e(:,1));
    qval = polyval(qp,e(:,2));
    rval = polyval(rp,e(:,3));

    dpval = polyval(dpp,e(:,1));
    dqval = polyval(dqp,e(:,2));
    drval = polyval(drp,e(:,3));
    
    % Normalization factor to ensure integration to one    
    fc = sqrt((2*pq(ii,1)+1)*2*(pq(ii,1)+pq(ii,2)+1)*2*(pq(ii,1)+pq(ii,2)+pq(ii,3)+2));

    fx(:,ii) = fc*(dpval.*qval.*rval.*de1(:,1));
    fy(:,ii) = fc*(dpval.*qval.*rval.*de1(:,2) + pval.*dqval.*rval.*de2(:,2));
    fz(:,ii) = fc*(dpval.*qval.*rval.*de1(:,3) + pval.*dqval.*rval.*de2(:,3) + pval.*qval.*drval);
end
fx = 2*fx;
fy = 2*fy;
fz = 2*fz;

function pq = pascalindex3d(p)
l=1;
for i=0:p
    for j=0:i
        for k=0:j
            pq(l,1)=i-j;
            pq(l,2)=j-k;
            pq(l,3)=k;
            l = l+1;
            if l>p 
               return;
            end
        end
    end
end
