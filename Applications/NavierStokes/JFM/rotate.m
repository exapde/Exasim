function R = rotate(n)
W = [0, -n(2), -n(3); n(2), 0, 0; n(3), 0, 0];
R = eye(3) + W + W*W/(1.0+n(1));
end

