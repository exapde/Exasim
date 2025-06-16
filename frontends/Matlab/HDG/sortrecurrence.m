function [b, e] = sortrecurrence(a)

n = length(a);
recurrence = zeros(1, n);

% Step 1: Count occurrences for each entry
for i = 1:n
    count = 0;
    for j = 1:n
        if a(i) == a(j)
            count = count + 1;
        end
    end
    recurrence(i) = count;
end

maxrec = 0;
for i = 1:n
    if recurrence(i) > maxrec
        maxrec = recurrence(i);
    end
end

b = zeros(1, n);
e = zeros(1, n);
c = zeros(1, n);
d = zeros(1, n);
m = 0;
for k = 1:maxrec
    j = 0;    
    for i = 1:n
        if recurrence(i) == k
            j = j + 1;
            c(j) = a(i);        
            d(j) = i;
        end
    end
    [c, ind] = simple_bubble_sort(c(1:j));
    d = d(ind);
    b(m+1:m+j) = c;
    e(m+1:m+j) = d;
    m = m + j;
    if m>=n
      break;
    end
end

if max(abs(b-a(e))) > 0
  error("something wrong");
end

end

function [b, ind] = simple_bubble_sort(a)
n = length(a);
b = a;
ind = 1:n;

for i = 1:n-1
    for j = 1:n-i
        if b(j) > b(j+1)
            % Swap values
            tmp = b(j);
            b(j) = b(j+1);
            b(j+1) = tmp;
            % Swap indices
            tmp_idx = ind(j);
            ind(j) = ind(j+1);
            ind(j+1) = tmp_idx;
        end
    end
end
end

