function ind = outerwall(p)

load midcurve.mat
n = size(mid,1);
d = (mid(:,1) - p(1)).^2 + (mid(:,2) - p(2)).^2;
[~,i] = min(d);
if i == 1
  [~, ~, position] = find_perpendicular_point(mid(i,1), mid(i,2), mid(i+1,1), mid(i+1,2), p(1), p(2));
elseif i==n
  [~, ~, position] = find_perpendicular_point(mid(i-1,1), mid(i-1,2), mid(i,1), mid(i,2), p(1), p(2));
else
  [~, ~, position1] = find_perpendicular_point(mid(i,1), mid(i,2), mid(i+1,1), mid(i+1,2), p(1), p(2));
  [~, ~, position2] = find_perpendicular_point(mid(i-1,1), mid(i-1,2), mid(i,1), mid(i,2), p(1), p(2));
  if (position1 == 1) || (position2==1)
    position = 1;
  else
    position = 0;
  end
end

if position==1
  ind = 1;
else 
  ind = 0;
end

