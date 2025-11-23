function d = parselogfile1(file)

d1 = parselogfile(file, 'hdgAssembleLinearSystem time:');
d2 = parselogfile(file, 'Matrix-vector product time:');
d3 = parselogfile(file, 'Applying preconditioner time:');
d4 = parselogfile(file, 'Orthgonalization time:');
d5 = parselogfile(file, 'Solution update time:');
d6 = parselogfile(file, 'GMRES time:');

d7 = parselogfile0(file, 'GMRES(', "Warning:");
d8 = parselogfile0(file, 'the tolerance', "Warning:");
d9 = parselogfile0(file, 'within', "Warning:");

% n = length(d2);
% d = [d1(1:n) d2 d3 d4 d5 d6(1:n) d7(1:n) d8(1:n) d9(1:n)];

d10 = parselogfile(file, 'Alpha:');
d11 = parselogfile(file, 'Original Norm:');
d12 = parselogfile(file, 'Updated Norm:');

A = [d10 d11 d12];
keep = true(size(A,1),1);  % mask for rows to keep
n = size(A,1);
for i = 2:n
    % If the first column value increases -> new group
    if A(i,1) > A(i-1,1)
        keep(i-1) = true;  % keep the previous row
    elseif A(i,1) == A(i-1,1)
        % same value -> keep both
        keep(i) = true;
    elseif A(i,1) < A(i-1,1)
        % decreasing trend: remove previous row(s) with larger first-column value
        keep(i-1) = false;
    end
end
B = A(keep,:);

n = length(d2);
[length(d1) length(d2) length(d3) length(d4) length(d5) length(d6) length(d7) length(d8) length(d9) length(d10) length(d11) length(d12)]
d = [d1(1:n) d2 d3 d4 d5 d6(1:n) d7(1:n) d8(1:n) d9(1:n) B(1:n,:)];

