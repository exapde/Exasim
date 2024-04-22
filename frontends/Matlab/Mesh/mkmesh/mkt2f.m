function [f,t2f,t2t] = mkt2f(t,elemtype)
%MKT2F Compute Face Connectivity and Triangle to Face Connetivity.
%   [F,T2F]=MKT2F(T)
%
%      T:         Triangle indices (NT,3)
%      ELEMTYPE:  Flag determining element type
%                 Flag = 0 tri/tet elements (default)
%                 Flag = 1 quad/hex elements
%
%      F:         Face connectivity (NF,4) (for boundary edges F(:,4)=0)
%      T2F:       Triangle to Face Connectivity (NT,3)
%
%   See also: MKT2T.

% npv : number of nodes per volume element
% nfv : number of faces per volume element
% npf : number of nodes per face element

if nargin<2, elemtype=0; end

[nt,npv]=size(t);
if npv==2 % 1D
    dim=1;
    nfv=2;
    npf=1;
else
    if elemtype==0 % tri/tet elements
        dim=size(t,2)-1;        
        nfv=dim+1;
        npf = dim;
    else % quad/hex elements
        dim=log2(size(t,2));        
        nfv=2*dim;
        npf=2*(dim-1);
    end
end

switch dim
    case 1
        face=[1;2];
    case 2
        if elemtype==0
            face=[[2,3];[3,1];[1,2]];
        elseif elemtype==1
            face=[[1,2];[2,3];[3,4];[4,1]];
        end
    case 3
        if elemtype==0
            face=[[2,3,4];[1,4,3];[1,2,4];[1,3,2]];
        elseif elemtype==1
            face=[[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]];
        end
    otherwise
        error('Only can handle dim=1, dim=2 or dim=3');
end


t2t = mkt2t(t,elemtype);
nb = sum(sum(t2t <= 0));
f = zeros((nfv*nt+nb)/2,npf+2);
t2f = zeros(nt,nfv);
jf = 0;
for i=1:nt
    for j=1:nfv
        if t2t(i,j) > i || t2t(i,j) <=0
            ie = t2t(i,j);
            jf = jf + 1;
            
            f(jf,1:npf) = t(i,face(j,:));
            f(jf,npf+1) = i;
            f(jf,npf+2) = ie;
            t2f(i,j) = jf;
            
            if ie > 0
                %k = sum(reshape(t(ie,face),[nfv npf]),2)-sum(f(jf,1:npf))==0;                                                            
                %t2f(ie,k) = jf;                
                a = sort(reshape(t(ie,face),[nfv npf]),2);
                b = sort(f(jf,1:npf));                
                k = sum(abs(a-repmat(b,[nfv 1])),2)==0;
                t2f(ie,k) = jf;                                
            end                        
        end
    end
end


% Reorder faces - First interior then boundary

% nf = size(f,1);
% [a,mp] = sort(f(:,npf+2) == 0);
% f = f(mp,:);
% 
% a = zeros(nf,1);
% a(mp) = (1:nf)';
% t2f = reshape(t2f,nfv*nt,1);
% ii = find(t2f > 0);
% t2f(ii) = a(t2f(ii));
% ii = find(t2f < 0);
% t2f(ii) = -a(-t2f(ii));
% t2f = reshape(t2f,nt,nfv);



