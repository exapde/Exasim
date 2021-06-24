function [facecon,f,t2f,t2t,t] = mkconnectivity(p, t, porder, elemtype, bndexpr, periodicexpr)

disp('run mkt2f...');           
[f,t2f,t2t] = mkt2f(t,elemtype);

disp('run setbndnbrs...');           
f = setbndnbrs(p,f,bndexpr);    

disp('Create facecon...');
dim = size(p,2);
[philocvl,philocfc,~,~,perm] = localbasis(porder,dim,elemtype);
philocfc = philocfc{1};
perm = cell2mat(perm);

nf = size(f,1);
np = max(perm(:));
npf = size(perm,1);

facecon = zeros(npf,2,nf);
for i = 1:nf
    fn = f(i,1:end-2);
    fi = f(i,end-1:end); % obtain two elements sharing the same face i      
    if fi(2)>0           % face i is an interior face       
        kf = t2f(fi,:);         % obtain neighboring faces 
        i1 = kf(1,:)==i;  % obtain the index of face i in the 1st element
        i2 = kf(2,:)==i;  % obtain the index of face i in the 2nd element                                            
        
        pf = philocfc*p(fn,:); % nodal points on face i    
        pf = snap(pf);
        t1 = t(fi(1),:);
        p1 = philocvl*p(t1,:);
        p1 = p1(perm(:,i1),:);
        p1 = snap(p1);     
        j1 = xiny(p1,pf);                
        t2 = t(fi(2),:);
        p2 = philocvl*p(t2,:);      
        p2 = p2(perm(:,i2),:);
        p2 = snap(p2);                                       
        j2 = xiny(p2,pf);                               
        
        facecon(:,1,i) = (fi(1)-1)*np+perm(j1,i1);
        facecon(:,2,i) = (fi(2)-1)*np+perm(j2,i2);        
    else % face i is a boundary face
        kf = t2f(fi(1),:); % obtain neighboring faces 
        i1 = kf(1,:)==i;  % obtain the index of face i in the 1st element   
        
        pf = philocfc*p(fn,:); % nodal points on face i    
        pf = snap(pf);
        t1 = t(fi(1),:);
        p1 = philocvl*p(t1,:);
        p1 = p1(perm(:,i1),:);
        p1 = snap(p1);     
        j1 = xiny(p1,pf);                
        
        facecon(:,1,i) = (fi(1)-1)*np+perm(j1,i1);        
        facecon(:,2,i) = facecon(:,1,i);
    end        
end 
%facecon = reshape(facecon,[npf*nf 2])'; 

if isempty(periodicexpr)==0    
    % get dgnodes
    npv = size(philocvl,1);
    [ne, nnv] = size(t);
    nd = size(p,2);
    dgnodes = zeros(npv,ne,nd);
    for d=1:nd
      for n=1:nnv
        dp=philocvl(:,n)*p(t(:,n),d)';
        dgnodes(:,:,d)=dgnodes(:,:,d)+dp;
      end
    end    
    dgnodes = reshape(dgnodes,[npv*ne,nd]);
    dgnodes = snap(dgnodes);
    
    disp('Handle periodic boundaries...');            
    [t,f,t2f,t2t,facecon] = periodic(p,t,f,t2f,t2t,facecon,periodicexpr);    
    
    % fix facecon so that DOFs on periodic faces are matched.
    nf = size(f,1);
    nperiodic = size(periodicexpr,1);
    for i = 1:nf
        x1 = dgnodes(facecon(:,1,i),:); % dgnodes on the 1st element
        x2 = dgnodes(facecon(:,2,i),:); % dgnodes on the 2nd element
        if max(abs(x1(:)-x2(:)))>1e-8  % periodic face
            for j = 1:nperiodic
                p = x1;
                y1 = eval(periodicexpr{j,2}); 
                p = x2;
                y2 = eval(periodicexpr{j,4}); 
                % sort y2 to match it to y1
                ind = xiny(y1, y2);
                if min(ind)>0
                    facecon(:,2,i) = facecon(ind,2,i);
                end                  
            end
        end
    end
     
    for i = 1:nf
        fi = f(i,end-1:end); % obtain two elements sharing the same face i           
        if fi(2)>0           % face i is an interior face     
%             [i fi]
%             [facecon(:,1,fi(1)) facecon(:,2,fi(1))]
%             [dgnodes(facecon(:,1,fi(1)),:) dgnodes(facecon(:,2,fi(1)),:)]
%             d = dgnodes(facecon(:,1,i),:)-dgnodes(facecon(:,2,i),:);
%             if max(abs(d(:)))>1e-10
%                 disp([i fi]);
%                 disp([dgnodes(facecon(:,1,i),:) dgnodes(facecon(:,2,i),:)]);
%                 error('periodic function is wrong');
%             end
            
            kf = t2f(fi,:);         % obtain neighboring faces 
            if ~ismember(i,kf(1,:)) % obtain the index of face i in the 1st element
                error('periodic function is wrong');
            end            
            if ~ismember(i,kf(2,:)) % obtain the index of face i in the 1st element
                error('periodic function is wrong');
            end                        
        else % boundary faces
            kf = t2f(fi(1),:);         % obtain neighboring faces 
            if ~ismember(i,kf(1,:)) % obtain the index of face i in the 1st element
                error('periodic function is wrong');
            end                        
        end    
    end
    
end

if size(t2f,1) ~= size(t,1)
    error('something is wrong');
end                        
if size(t2t,1) ~= size(t,1)
    error('something is wrong');
end        
if size(facecon,3) ~= size(f,1)
    error('something is wrong');
end                        
            
function f=setbndnbrs(p0,f,bndexpr)
%SETBNDNBRS Set Boundary Marker for the Boundary Faces.
%   F=SETBNDNBRS(P0,F,BNDEXPR)
%
%      P0:        Node positions (NP,2)
%      F:         Face Array (NF,4)
%      BNDEXPR:   Cell Array of boundary expressions. The 
%                 number of elements in BNDEXPR determines 
%                 the number of different boundaries
%
%   Example: (Setting boundary types for a unit square mesh - 4 types)
%      bndexpr = {'all(p(:,2)<1e-3)','all(p(:,1)>1-1e-3)', ...
%                 'all(p(:,2)>1-1e-3)','all(p(:,1)<1e-3)'};     
%      f = setbndnbrs(p,f,bndexpr);
%
%   Example: (Setting boundary types for the unit circle - 1 type)
%      bndexpr = bndexpr = {'all(sqrt(sum(p.^2,2))>1-1e-3)'}; 
%      f = setbndnbrs(p,f,bndexpr);
%
dim = size(p0,2);

%[i,foo]=find(f==0);
i=find(f(:,end)==0);

for ii=i'
  p=p0(f(ii,1:dim),:);
  
  found=false;
  for jj=1:length(bndexpr)
    if eval(bndexpr{jj})
      found=true;
      bnd=jj;
      break;
    end
  end  
  
  if ~found
    error('Strange boundary.');
  end
  
  f(ii,end)=-bnd;
end


