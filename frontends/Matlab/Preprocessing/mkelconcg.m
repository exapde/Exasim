% function [cgnodes, cgelcon, ncg] = mkelconcg(dgnodes)
% % Inputs:
% %   dgnodes : [npe, dim, ne] array of discontinuous Galerkin node coordinates
% % Outputs:
% %   cgnodes : [max_nodes, dim] array of unique continuous Galerkin node coordinates
% %   cgelcon : [npe, ne] connectivity array mapping local DG node to CG node index
% %   ncg     : number of unique CG nodes
% 
% tol = 1e-8;
% [npe, dim, ne] = size(dgnodes);
% 
% % Preallocate with max possible number of nodes (npe * ne)
% max_nodes = npe * ne;
% cgnodes = zeros(max_nodes, dim);
% cgelcon = zeros(npe, ne);
% 
% ncg = 0;  % current number of unique CG nodes
% 
% for e = 1:ne
%     for a = 1:npe
%         node = reshape(dgnodes(a, :, e),[1 dim]);  % 1×dim
% 
% %         found = -1;
% %         for j = 1:ncg
% %             if max(abs(cgnodes(j, :) - node)) < tol
% %                 found = j;
% %                 break;
% %             end
% %         end
% 
%         if (ncg >= 1) 
%           diffs = abs(cgnodes(1:ncg,:) - node);                % compute absolute differences
%           dists = max(diffs, [], 2);                  % ?-norm across rows
%           found = find(dists < tol, 1);             % find first match
%           if isempty(found)
%               found = -1;
%           end        
%         else
%           found = -1; 
%         end
% 
% 
%         if found > 0
%             cgelcon(a, e) = found;
%         else
%             ncg = ncg + 1;
%             cgnodes(ncg, :) = node;
%             cgelcon(a, e) = ncg;
%         end
%     end
% end
% 
% % Trim cgnodes to the number of unique nodes
% cgnodes = cgnodes(1:ncg, :);
% end

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





