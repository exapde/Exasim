function mesh = joinBoundaries(mesh,b1,b2,ic)
%JOINBOUNDARIES Function to merge two boundaries into a periodic boundary,
%               by modifying only mesh.t2f and mesh.elcon.
%               NOTE: run this BEFORE calling preprocess.m.
%   mesh = joinBoundaries2(mesh,b1,b2,ic)
%
%     b1 = boundary 1
%     b2 = boundary 2 (periodic pair to boundary 1)
%     ic = which coordinate to use to identify paired faces (e.g. 2 if x-periodic)
%
%--------------------------------------------------------------------------
% REVISION HISTORY:
% When     Who               What
% 17Oct12  Hemant Chaurasia  Modified from joinBoundaries_old
%--------------------------------------------------------------------------

% Find matching faces to merge
[fb1 foo] = find(mesh.f(:,4)==-b1);
[fb2 foo] = find(mesh.f(:,4)==-b2);
if length(fb1)~=length(fb2); error('Can''t match faces.'); end
fmerge = zeros(length(fb1),2);
fc = (mesh.p(mesh.f(:,1),:)+mesh.p(mesh.f(:,2),:))/2; % centroid of all points
for i1=1:length(fb1)
    f1 = fb1(i1);
    fmerge(i1,1) = f1;

    % Find corresponding face on boundary b2
    fc1 = fc(f1,:); % face 1 centroidal point
    [ii foo] = find(abs(fc(fb2,ic)-fc1(ic))<2e-15);
    if length(ii)~=1; 
        error(['Matching error, i1 = ', num2str(i1)]); 
    end
    fmerge(i1,2) = fb2(ii);
end

% Delete paired b2 faces
newf = mesh.f;
newt2f = mesh.t2f;
for i=1:size(fmerge,1)
    newf(fmerge(i,1),4) = newf(fmerge(i,2),3);  % relabel f(:,4) using f(:,3) from fmerge(i,2)
    
    for j=1:size(newt2f,1)
        for k=1:size(newt2f,2);
            if newt2f(j,k)==fmerge(i,2)
                newt2f(j,k) = fmerge(i,1);      % relabel fmerge(i,2) --> fmerge(i,1) in t2f
            end
            if newt2f(j,k)>fmerge(i,2)
                newt2f(j,k) = newt2f(j,k)-1;    % decrement all face numbers > fmerge(i,2)
            end
        end
    end 
    
    % Delete face fmerge(i,2)
    newf = newf([1:fmerge(i,2)-1 fmerge(i,2)+1:end],:);
    for j=1:size(fmerge,1)
        for k=1:size(fmerge,2)
            if fmerge(j,k)>fmerge(i,2)
                fmerge(j,k) = fmerge(j,k)-1;  % update-decrement fmerge face numbers too
            end
        end
    end   
end

% Modify mesh.elcon

mesh.f = newf;
mesh.t2f = newt2f;
mesh.fcurved = zeros(length(newf),1);
mesh.nf = length(newf);