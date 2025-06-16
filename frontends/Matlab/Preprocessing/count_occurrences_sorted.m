function [b, c] = count_occurrences_sorted(a)
% COUNT_OCCURRENCES_SORTED  Count runs in a sorted vector and their positions
%   [b,c] = COUNT_OCCURRENCES_SORTED(a) returns two vectors of same size as
%   sorted input a:
%     b(i) = number of times a(i) appears in a
%     c(i) = the 1-based position of a(i) within its run of identical values
%
% Example:
%   a = [2 2 2 5 5 7 9 9];
%   [b,c] = count_occurrences_sorted(a);
%   % b = [3 3 3 2 2 1 2 2]
%   % c = [1 2 3 1 2 1 1 2]

n = numel(a);
b = zeros(size(a));
c = zeros(size(a));

i = 1;
while i <= n
    % find end of this run
    j = i + 1;
    while j <= n && a(j) == a(i)
        j = j + 1;
    end
    cnt = j - i;
    % fill b and c for positions i:(j-1)
    b(i:j-1) = cnt;
    c(i:j-1) = 1:cnt;
    % advance
    i = j;
end

end
