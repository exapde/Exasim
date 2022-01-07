function [fcolor,ncolor] = colorface(f,t2f,perm,faceind)
% break the face block faceind into smaller groups


nfe = size(perm,2);
a = ones(nfe,nfe);
b = [];
for i = 1:nfe
    for j = 1:nfe
        if j~=i
            if isempty(intersect(perm(:,i),perm(:,j)))
                % local face i and local face j are not connected
                a(i,j) = 0;
            end
        end
    end
    b = [b; setdiff(find(a(i,:)==1),i)];
end

%faceind = if1:1:if2;
%nf = if2-if1+1;
nf = length(faceind);
nb = size(b,2);
fb = zeros(nf,2*nb);
for i = 1:nf    
    fi = faceind(i);          % global face fi
    ei = f(fi,end-1:end); % two elements share face fi
    f1 = t2f(ei(1),:); % list of neighboring faces on element ei(1)
    i1 = f1==fi; % local index of face fi on element ei(1)
    b1 = b(i1,:); % local indices of faces connected to face fi
    fb(i,1:nb) = t2f(ei(1),b1); % list of global faces connected to global face fi on element ei(1)   
    if ei(2)>0
        f2 = t2f(ei(2),:); % list of neighboring faces on element ei(2)
        i2 = f2==fi; % local index of face fi on element ei(2)
        b2 = b(i2,:); % local indices of faces connected to face fi
        fb(i,nb+1:end) = t2f(ei(2),b2);  % list of global faces connected to global face fi on element ei(2)         
    end        
end

c = 1;
fcolor = zeros(1,nf);
ncolor = [0 0];
while (1)        
    if c == 1 % color number 1
        [color,inda] = setcolor(fb,faceind); % find all faces have the same color number 1
        ncolor(c+1) = ncolor(c) + length(color); % set the number of faces in the list
        fcolor(ncolor(c)+1:ncolor(c+1)) = color; % add those faces to the list
        ind = inda;
    else % color number c
        ind = setdiff(1:nf, inda); % remove all the previous faces        
        if ~isempty(ind)
            [color,in] = setcolor(fb(ind,:),faceind(ind));  % find all faces have the same color number c                        
            ncolor(c+1) = ncolor(c) + length(color); % set the number of faces in the list
            fcolor(ncolor(c)+1:ncolor(c+1)) = color; % add those faces to the list
            inda = [inda ind(in)]; % update the index list
        end
    end    
    if isempty(ind)
        break;
    else
        c = c + 1;
    end
end

check = 1;
if check==1
    npe = max(perm(:));
    npf = size(perm,1);
    for k = 1:c-1
        fc = fcolor(ncolor(k)+1:ncolor(k+1));        
        u1 = zeros(npf,length(fc));
        u2 = zeros(npf,length(fc));
        for i=1:length(fc)
            fi = fc(i);           % global face fi
            ei = f(fi,end-1:end); % two elements share face fi
            f1 = t2f(ei(1),:); % list of neighboring faces on element ei(1)
            i1 = f1==fi; % local index of face fi on element ei(1)
            u1(:,i) = (ei(1)-1)*npe+perm(:,i1); % DOFs of udg on element ei(1)
            if ei(2)>0
                f2 = t2f(ei(2),:); % list of neighboring faces on element ei(2)
                i2 = f2==fi; % local index of face fi on element ei(2)                
                u2(:,i) = (ei(2)-1)*npe+perm(:,i2); % DOFs of udg on element ei(1)
            end                                
        end                  
        if length(unique(u1(:))) ~= numel(u1)                              
            error('something wrong');
        end
        u2 = u2(u2>0);        
        if length(unique(u2(:))) ~= numel(u2)                              
            error('something wrong');
        end        
    end
end


function [flist,ind] = setcolor(fb,faceind)

nf = size(fb,1);
flist = ones(1,nf); % list of the faces having the same color
flist(1) = faceind(1); % starting with the first face
ind = ones(1,nf);
k = 1;
for i = 2:nf % for each face i           
    if isempty(intersect(fb(i,:),flist(1:k)))
        % if any neighboring faces of face i are not in the list
        k = k + 1;        
        flist(k) = faceind(i); % add face i to the list       
        ind(k) = i;
    end    
end
flist = flist(1:k);
ind = ind(1:k);


