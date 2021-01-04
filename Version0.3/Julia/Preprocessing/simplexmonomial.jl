function monomial1d(x::Array{Float64,2},p::Int)::Tuple{Array{Float64,2},Array{Float64,2}}

n = length(x);
v0 = zeros(Float64,n,1);
v1 = ones(Float64,n,1);

if (p==0)
    f = v1;
    fx = v0;
elseif (p==1)
    f = [v1 x];
    fx = [v0 v1];
elseif (p==2)
    f = [v1 x x.^2];
    fx = [v0 v1 2*x];
elseif (p==3)
    f = [v1 x x.^2 x.^3];
    fx = [v0 v1 2*x 3*x.^2];
elseif (p==4)
    f = [v1 x x.^2 x.^3 x.^4];
    fx = [v0 v1 2*x 3*x.^2 4*x.^3];
else
    error("polynomial degree must be less than 5");
end

return f, fx;

end

function monomial2d(xy::Array{Float64,2},p::Int)::Tuple{Array{Float64,2},Array{Float64,2},Array{Float64,2}}

n = size(xy,1);
v0 = zeros(Float64,n,1);
v1 = ones(Float64,n,1);

x = xy[:,1];
y = xy[:,2];

if (p==0)
    f = v1;
    fx = v0;
    fy = v0;
elseif (p==1)
    f = [v1 x y];
    fx = [v0 v1 v0];
    fy = [v0 v0 v1];
elseif (p==2)
    f = [v1 x y x.*y x.^2 y.^2];
    fx = [v0 v1 v0 y 2*x v0];
    fy = [v0 v0 v1 x v0 2*y];
elseif (p==3)
    f = [v1 x y x.*y x.^2 y.^2 x.*y.^2 (x.^2).*y x.^3 y.^3];
    fx = [v0 v1 v0 y 2*x v0 y.^2 2*x.*y 3*x.^2 v0];
    fy = [v0 v0 v1 x v0 2*y 2*x.*y x.^2 v0 3*y.^2];
elseif (p==4)
    f = [v1 x y x.*y x.^2 y.^2 x.*y.^2 (x.^2).*y x.^3 y.^3 x.*y.^3 (x.^2).*y.^2 (x.^3).*y x.^4 y.^4];
    fx = [v0 v1 v0 y 2*x v0 y.^2 2*x.*y 3*x.^2 v0 y.^3 2*x.*y.^2 3*(x.^2).*y 4*x.^3 v0];
    fy = [v0 v0 v1 x v0 2*y 2*x.*y x.^2 v0 3*y.^2 3*x.*y.^2 2*(x.^2).*y x.^3 v0 4*y.^3];
else
    error("polynomial degree must be less than 5");
end

return f, fx, fy;

end

function monomial3d(xyz::Array{Float64,2},p::Int)::Tuple{Array{Float64,2},Array{Float64,2},Array{Float64,2},Array{Float64,2}}

n = size(xyz,1);
v0 = zeros(Float64,n,1);
v1 = ones(Float64,n,1);

x = xyz[:,1];
y = xyz[:,2];
z = xyz[:,3];

if (p==0)
    f = v1;
    fx = v0;
    fy = v0;
    fz = v0;
elseif (p==1)
    f = [v1 x y z];
    fx = [v0 v1 v0 v0];
    fy = [v0 v0 v1 v0];
    fz = [v0 v0 v0 v1];
elseif (p==2)
    f = [v1 x y z x.*y x.*z y.*z x.^2 y.^2 z.^2];
    fx = [v0 v1 v0 v0 y z v0 2*x v0 v0];
    fy = [v0 v0 v1 v0 x v0 z v0 2*y v0];
    fz = [v0 v0 v0 v1 v0 x y v0 v0 2*z];
elseif (p==3)
    f = [v1 x y z x.*y x.*z y.*z x.^2 y.^2 z.^2 x.*y.*z x.*y.^2 x.*z.^2 (x.^2).*y (x.^2).*z y.*z.^2 (y.^2).*z x.^3 y.^3 z.^3];
    fx = [v0 v1 v0 v0 y z v0 2*x v0 v0 y.*z y.^2 z.^2 2*x.*y 2*x.*z v0 v0 3*x.^2 v0 v0];
    fy = [v0 v0 v1 v0 x v0 z v0 2*y v0 x.*z 2*x.*y v0 x.^2 v0 z.^2 2*y.*z v0 3*y.^2 v0];
    fz = [v0 v0 v0 v1 v0 x y v0 v0 2*z x.*y v0 2*x.*z v0 x.^2 2*y.*z y.^2 v0 v0 3*z.^2];
elseif (p==4)
    f = [v1 x y z x.*y x.*z y.*z x.^2 y.^2 z.^2 ...
         x.*y.*z x.*y.^2 x.*z.^2 (x.^2).*y (x.^2).*z y.*z.^2 (y.^2).*z x.^3 y.^3 z.^3 ...
         x.*y.*z.^2 (x.*y.^2).*z (x.^2).*y.*z (x.^2).*y.^2 (x.^2).*z.^2 (y.^2).*z.^2 ...
         (x.^3).*y (x.^3).*z x.*y.^3 (y.^3).*z x.*z.^3 y.*z.^3 x.^4 y.^4 z.^4];
    fx = [v0 v1 v0 v0 y z v0 2*x v0 v0 y.*z y.^2 z.^2 2*x.*y 2*x.*z v0 v0 3*x.^2 v0 v0...
          y.*z.^2 (y.^2).*z 2*x.*y.*z 2*x.*y.^2 2*x.*z.^2 v0 ...
          3*(x.^2).*y 3*(x.^2).*z y.^3 v0 z.^3 v0 4*x.^3 v0 v0];
    fy = [v0 v0 v1 v0 x v0 z v0 2*y v0 x.*z 2*x.*y v0 x.^2 v0 z.^2 2*y.*z v0 3*y.^2 v0...
          x.*z.^2 2*x.*y.*z (x.^2).*z 2*(x.^2).*y v0 2*y.*z.^2 ...
          x.^3 v0 3*x.*y.^2 3*(y.^2).*z v0 z.^3 v0 4*y.^3 v0];
    fz = [v0 v0 v0 v1 v0 x y v0 v0 2*z x.*y v0 2*x.*z v0 x.^2 2*y.*z y.^2 v0 v0 3*z.^2 ...
          2*x.*y.*z x.*y.^2 (x.^2).*y v0 2*(x.^2).*z 2*(y.^2).*z ...
          v0 x.^3 v0 y.^3 3*x.*z.^2 3*y.*z.^2 v0 v0 4*z.^3];
else
    error("polynomial degree must be less than 5");
end

return f, fx, fy, fz;

end

function simplexmonomial(x::Array{Float64,2},porder::Int)::Tuple{Array{Float64,2},Array{Float64,2},Array{Float64,2},Array{Float64,2}}

dim=size(x,2);

if dim==1
     f,fx=monomial1d(x,porder);      # 1D
     fy=[0.0 0.0];
     fz=[0.0 0.0];
elseif dim==2
     f,fx,fy=monomial2d(x,porder);   # 2D
     fz=[0.0 0.0];
elseif dim==3
     f,fx,fy,fz=monomial3d(x,porder);# 3D
else
     error("Only can handle dim=1, dim=2 or dim=3");
end

return f, fx, fy, fz;

end
