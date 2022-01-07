function cubemesh(m,n,o,elemtype)
# Generate mesh for unit cube

p = [[x; y; z] for x = 0.0:1.0/Float64(m):1.0, y = 0.0:1.0/Float64(n):1.0, z = 0.0:1.0/Float64(o):1.0][:]
p = hcat(p...);

m = m+1;
n = n+1;
o = o+1;

if elemtype==0
    t0=[[5 6 7 3]; [5 6 2 3]; [5 1 2 3]; [6 8 7 4]; [6 7 4 3]; [6 2 4 3]];

    map=[1 2];
    map=[map  map.+m];
    map=[map map.+m*n];
    t=map[t0];
    t=kron(t,ones(m-1,1))+kron(ones(size(t)),collect(0:m-2));
    t=kron(t,ones(n-1,1))+kron(ones(size(t)),collect(0:n-2).*m);
    t=kron(t,ones(o-1,1))+kron(ones(size(t)),collect(0:o-2).*(m*n));
else
    t = [1 2 m+2 m+1 m*n+1 m*n+2 m*n+m+2 m*n+m+1];
    t=kron(t,ones(o-1,1))+kron(ones(size(t)),collect(0:o-2).*(m*n));
    t=kron(t,ones(n-1,1))+kron(ones(size(t)),collect(0:n-2).*m);
    t=kron(t,ones(m-1,1))+kron(ones(size(t)),collect(0:m-2));
end

t = Int.(t');

return p, t

end
