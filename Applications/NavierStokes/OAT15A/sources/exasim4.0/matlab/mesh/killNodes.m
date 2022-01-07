
function [p2,t2] = killNodes(p,t,ih)

% Function [p2,t2] = killNodes(p,t,ih)
% Removes specified nodes' associated triangles, then remeshes the 
% resulting void using Delaunay.

p2=p; t2=t;
for i=1:length(ih)
    adjacentTriangles = find(sum(ismember(t2,ih(i)),2)>0);  % triangles containing this hanging node
    numAdjacentTriangles = length(adjacentTriangles);
    
    otherVerticesInTriangle = zeros(length(adjacentTriangles),2);
    for ii=1:numAdjacentTriangles
        otherVerticesInTriangle(ii,:) = setdiff(t2(adjacentTriangles(ii),:),ih(i)); % [numAdjacentTriangles,2]
    end
    verticesForTriangulation = reshape(otherVerticesInTriangle',1,[]);   % [2*numAdjacentTriangles]
    [uniqueVerticesForTriangulation,map1,map2] = unique(verticesForTriangulation);
    numUniqueVerticesForTriangulation = length(uniqueVerticesForTriangulation);
    
    constraintsForDelaunay = reshape(map2,2,numAdjacentTriangles)';
    
    if numAdjacentTriangles < numUniqueVerticesForTriangulation
        % The only way this can happen with simplex 2d elements is if the
        % hanging node is on the boundary. Add final constraint for
        % boundary edge.
        
        j = 1;
        for ii=1:numUniqueVerticesForTriangulation
            if sum(sum(constraintsForDelaunay == ii)) < 2
                % This vertex has been used only once for the Delaunay
                % constraints. it should be used again for the boundary
                % edge.
                if j > 2
                    error('It looks like there are at least three vertices that have been used only once for the Delaunay constraints.')
                end
                constraintsForDelaunay(numAdjacentTriangles+1,j) = ii;
                j = j+1;
            end
        end
    end
    
    DT = delaunayTriangulation(p(uniqueVerticesForTriangulation,:),constraintsForDelaunay);
    
    interiorTriangles = isInterior(DT);
    
    tvoid = DT.ConnectivityList(interiorTriangles,:);
    
    tvoid = uniqueVerticesForTriangulation(tvoid);
    
    if size(tvoid,2)==3                      % add new triangles
        t2 = [t2; tvoid];                
    else
        t2 = [t2; tvoid'];           
    end
    t2 = t2(setdiff(1:size(t2,1),adjacentTriangles),:);     % remove old triangles
end

end
