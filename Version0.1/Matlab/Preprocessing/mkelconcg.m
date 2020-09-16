function [cgnodes,cgelcon] = mkelconcg(dgnodes)

% remove duplicate nodes in mesh.p1
dim = size(dgnodes,2);
[ns,dim,nt] = size(dgnodes(:,1:dim,:));
A = reshape(permute(dgnodes(:,1:dim,:),[1,3,2]),[ns*nt,dim]);

% B = flip(A,1);
% [~,I]=unique(B,'rows'); 
% B = B(sort(I),:);
% B = flip(B,1);
% [~,b]=ismember(A,B,'rows');
[B,~,b] = unique(A,'rows'); 
%[~,b]=ismember(A,B,'rows');

% CG mesh
cgnodes = B;
cgelcon = reshape(b,[ns nt]);





