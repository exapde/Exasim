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
