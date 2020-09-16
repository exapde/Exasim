function ent2elem = mkcgent2elem(dgnodes, cgnodes, cgelcon, colent2elem, rowent2elem)

% spatial dimension
nd = size(dgnodes,2);

% number of elements connected to CG nodes
ent2nelem = rowent2elem(2:end)-rowent2elem(1:end-1); 

% number of CG nodes
nent = length(ent2nelem);

% maximum number of elements
maxnelem = max(ent2nelem); 

% storing list of elements for each CG node
ent2elem = zeros(maxnelem,nent);

for i = 1:nent % for CG node i
    elem = colent2elem((rowent2elem(i)+1):rowent2elem(i+1)); % elements connected to CG node i
    nelem = ent2nelem(i);  % number of elements connected to CG node i      
    if nelem<maxnelem                           
        % add neighboring elements until more than or equal to maxnelem
        nbelem = elem2nbelem(cgelcon, colent2elem, rowent2elem, elem);
        while (length(nbelem) < maxnelem)              
            nbelem = elem2nbelem(cgelcon, colent2elem, rowent2elem, nbelem);
        end            
        
        % make ent2elem(:,i) for CG node i
        if (length(nbelem) == maxnelem)        
            ent2elem(:,i) = nbelem;
        else
            % remove some elements     
            elemd = setdiff(nbelem,elem);
            dgx = reshape(mean(dgnodes(:,:,elemd),1),nd,[])';            
            dix = dgx-cgnodes(i,:);
            dix = sum(dix.^2,2);
            [~,ind] = sort(dix);
            ent2elem(:,i) = [elem; elemd(ind(1:(maxnelem-numel(elem))))];
        end
    else
        ent2elem(:,i) = elem;
    end    
end

