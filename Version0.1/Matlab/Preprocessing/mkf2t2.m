function [f2t, t2t] = mkf2t2(t,elemtype,nd)

[~,ne] = size(t);

% get local faces of an element
face = getelemface(nd,elemtype);
[nvf,nfe] = size(face);

% sort faces for all elements to make two identical faces have the same numbering order
tf = sort(reshape(t(face,:),[nvf nfe*ne]),1);

[tf, jx] = sortrows(tf'); % do this so that two identical faces are next each other
tf = tf';
dx = tf(:,2:end)-tf(:,1:end-1); % do this to find two identical faces
dx = sum(dx.^2,1);           % identical faces corresponds to zero entries
in1 = find(dx==0);           % interior global faces connected to 1st elements (1st identical faces)
in2 = in1+1;                 % interior global faces connected to 2nd elements (2nd identical faces)
in0 = setdiff((1:length(jx))',unique([in1(:); in2(:)])); % boundary global faces connected to 1st elements

nf = length(in0)+length(in1);
f2t = zeros(4,nf);   % allocate memory 
t2t = zeros(nfe,ne); % allocate memory 

% interior faces
e1 = ceil(jx(in1)/nfe);      % 1st elements 
l1 = jx(in1) - (e1-1)*nfe;   % 1st local faces
e2 = ceil(jx(in2)/nfe);      % 2nd elements 
l2 = jx(in2) - (e2-1)*nfe;   % 2nd local faces
g = 1:length(in1);           % indices for interior faces
f2t(1,g) = e1;
f2t(2,g) = l1;
f2t(3,g) = e2;
f2t(4,g) = l2;                  
for i = 1:length(e1)
    t2t(l1(i),e1(i)) = e2(i);
    t2t(l2(i),e2(i)) = e1(i);
end

% boundary faces
e1 = ceil(jx(in0)/nfe);      % 1st element 
l1 = jx(in0) - (e1-1)*nfe;   % 1st local face
g = (length(in1)+1):nf;      % indices for interior faces
f2t(1,g) = e1;
f2t(2,g) = l1;


% % Python: c, ia, ic= np.unique(a,return_index=True,return_inverse=True, axis=0)
% 
% % do this to find two elements sharing the same face
% [~,ix,jx] = unique(tf','rows'); % 
% [ax, bx] = sort(jx); % do this so that two neigboring elements are next each other
%                      % ax -> indices of the global faces
%                      % bx -> indices of the elements and local faces
%                      
% na = length(ax);  % should be equal nfe*ne
% nf = length(ix);  % number of global unique faces
% f2t = zeros(4,nf); % allocate memory 
% t2t = zeros(nfe,ne);
% 
% dx = diff(ax);
% in1 = find(dx==0);           % interior global faces connected to 1st elements
% in2 = in1+1;                 % interior global faces connected to 2nd elements
% e1 = ceil(bx(in1)/nfe);      % 1st elements 
% l1 = bx(in1) - (e1-1)*nfe;   % 1st local faces
% e2 = ceil(bx(in2)/nfe);      % 2nd elements 
% l2 = bx(in2) - (e2-1)*nfe;   % 2nd local faces
% g = ax(in1);
% f2t(1,g) = e1;
% f2t(2,g) = l1;
% f2t(3,g) = e2;
% f2t(4,g) = l2;                  
% for i = 1:length(e1)
%     t2t(l1(i),e1(i)) = e2(i);
%     t2t(l2(i),e2(i)) = e1(i);
% end
% 
% in1 = setdiff((1:na)',unique([in1(:); in2(:)])); % boundary global faces connected to 1st elements
% e1 = ceil(bx(in1)/nfe);      % 1st element 
% l1 = bx(in1) - (e1-1)*nfe;   % 1st local face
% g = ax(in1);
% f2t(1,g) = e1;
% f2t(2,g) = l1;

% k = 1;
% while k<=na    
%     g = ax(k); % global face g
%     
%     if k==na % handle the special case to avoid checking if ax(k)==ax(k+1)
%         e1 = ceil(bx(k)/nfe);
%         l1 = bx(k) - (e1-1)*nfe;
%         f2t(1,g) = e1;
%         f2t(2,g) = l1;        
%         break;
%     end
%         
%     if ax(k)==ax(k+1) % this global face is connected to two elements (interior face)      
%         e1 = ceil(bx(k)/nfe);      % 1st element 
%         l1 = bx(k) - (e1-1)*nfe;   % 1st local face
%         e2 = ceil(bx(k+1)/nfe);    % 2nd element 
%         l2 = bx(k+1) - (e2-1)*nfe; % 2nd local face
%         f2t(1,g) = e1;
%         f2t(2,g) = l1;
%         f2t(3,g) = e2;
%         f2t(4,g) = l2;                  
%         t2t(l1,e1) = e2;
%         t2t(l2,e2) = e1;
%         k = k + 2; % increase by 2, because there are two elements connected to this face     
%     else  % this global face is connected to one element (boundary face)              
%         e1 = ceil(bx(k)/nfe);    % 1st element 
%         l1 = bx(k) - (e1-1)*nfe; % 1st local face
%         f2t(1,g) = e1;
%         f2t(2,g) = l1;
%         k = k + 1; % inceease by 1, because there is only one element connected to this face
%     end
% end

% tf = reshape(tf,[nvf nfe ne]);
% for i = 1:nf
%     if f2t(3,i)>0
%         e = tf(:,f2t(2,i),f2t(1,i))-tf(:,f2t(4,i),f2t(3,i));
%         if max(abs(e))>0
%             error("something wrong");
%         end
%     end
% end
