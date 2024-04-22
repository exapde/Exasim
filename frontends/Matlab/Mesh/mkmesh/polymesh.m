function [p,t]=polymesh(pvs,closed,regions,hpars,hfunc,varargin)
%POLYMESH Generate triangular mesh of polygon region.
%   [P,T]=POLYMESH(PVS,CLOSED,REGIONS,HPARS,HFUNC,HFUNCPARS)
%
%   Example: Square with hole
%     n=32;
%     phi=2*pi*(1:n)'/n;
%     pv1=[cos(phi),sin(phi)];
%     pv2=[-2,-2; 2,-2; 2,2; -2,2];
%     [p,t]=polymesh({pv1,pv2},[1,1],[0,1;1,0]);
%     simpplot(p,t)

if nargin<4, hpars=[1e3,1.3]; end

gtot=[];
for ii=1:length(pvs)
  pv=pvs{ii};
  np=size(pv,1);
  if closed(ii)
    g=[2*ones(1,np);pv(:,1)';pv([2:end,1],1)';pv(:,2)';pv([2:end,1],2)';
       regions(ii,1)*ones(1,np);regions(ii,2)*ones(1,np)];
  else
    g=[2*ones(1,np-1);pv([1:end-1],1)';pv([2:end],1)';pv([1:end-1],2)';pv([2:end],2)';
       regions(ii,1)*ones(1,np-1);regions(ii,2)*ones(1,np-1)];
  end
  gtot=[gtot,g];
end

[p,e,t]=initmesh(gtot,'hgrad',hpars(2),'hmax',hpars(1));

if nargin>=5
  while 1
    area=pdetrg(p,t);
    pmid=(p(:,t(1,:))+p(:,t(2,:))+p(:,t(3,:)))/3;
    hmid=feval(hfunc,pmid',varargin{:});
    it=find(area'>2*hmid.^2*sqrt(3)/4);
    if isempty(it)
      break;
    end
    [p,e,t]=refinemesh(gtot,p,e,t,[],it(:),'regular');
  end
  p=jigglemesh(p,e,t,'opt','mean','iter',100);
end

p=p';
t=t(1:3,:)';
