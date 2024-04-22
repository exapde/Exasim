function t2t = mkt2t(t,elemtype)
%MKT2T Compute Element to Element Connectivity.
%   T2T=MKT2T(T,MESHTYPE)
%
%      T:         Triangle indices (NT,3)
%      ELEMTYPE:  Flag determining element type
%                 Flag = 0 tri/tet elements (default)
%                 Flag = 1 quad/hex elements
%
%      T2T:       Triangle to Trangle Connectivity (NT,3)
%                 T2T(IT,IN) is the trangle that shares an edge
%                 with triangle IT and does nont contain node T(IT,IN).
%                 When an element is adjacent to the boundary the
%                 corresponding entry in T2T is set to zero
%

% npv : number of nodes per volume element
% nfv : number of faces per volume element
% npf : number of nodes per face element

if nargin<2, elemtype=0; end

[nt,npv]=size(t);

if npv==2
    dim=1;
    nfv=2;
else
    if elemtype==0 % tri/tet elements
        dim=size(t,2)-1;    
        nfv=dim+1;        
    elseif elemtype==1 % quad/hex elements
        dim=log2(size(t,2));   
        nfv=2*dim;
    end
end

switch dim
    case 1
        faces=[t(:,1)
               t(:,2)];
    case 2
        if elemtype==0
            faces=[t(:,[2,3])
                   t(:,[3,1])
                   t(:,[1,2])];
        elseif elemtype==1
            faces=[t(:,[1,2])
                   t(:,[2,3])
                   t(:,[3,4])
                   t(:,[4,1])];
        end
    case 3
        if elemtype==0
            faces=[t(:,[2,3,4])
                   t(:,[1,4,3])
                   t(:,[1,2,4])
                   t(:,[1,3,2])];
        elseif elemtype==1
            faces=[t(:,[1,4,3,2])
                   t(:,[5,6,7,8])
                   t(:,[1,2,6,5])
                   t(:,[3,4,8,7])
                   t(:,[2,3,7,6])
                   t(:,[4,1,5,8])];
        end
    otherwise
        error('Only can handle dim=1, dim=2 or dim=3');
end

ts=[repmat(int32(1:nt),1,npv); kron(int32(1:npv),ones(1,nt,'int32'))]';

faces=sort(faces,2);
%[faces(1:10,:); faces(end-10:end,:)]
[foa,fob,jx]=unique(faces,'rows');
[jx,ix]=sort(jx);
ts=ts(ix,:);
ix=find(diff(jx)==0);
ts1=ts(ix,:);
ts2=ts(ix+1,:);

t2t=zeros(nt,nfv,'int32');
t2t(ts1(:,1)+nt*(ts1(:,2)-1))=ts2(:,1);
t2t(ts2(:,1)+nt*(ts2(:,2)-1))=ts1(:,1);


