function legendrepolynomial(x,p)

n = length(x);
v0 = zeros(Float64,n,1);
v1 = ones(Float64,n,1);

x = 2.0*reshape(x,n,1) .- 1.0;

if (p==0)
    f = v1;
    fx = v0;
elseif (p==1)
    f = [v1 x];
    fx = [v0 v1];
elseif (p==2)
    f = [v1 x 0.5*(3*(x.^2) .- 1)];
    fx = [v0 v1 3*x];
elseif (p==3)
    f = [v1 x 0.5*(3*(x.^2) .- 1) 0.5*(5*(x.^3) .- 3*x)];
    fx = [v0 v1 3*x 7.5*(x.^2) .- 1.5];
elseif (p==4)
    f = [v1 x 0.5*(3*(x.^2) .- 1) 0.5*(5*(x.^3) .- 3*x) (1.0/8.0)*(35*(x.^4) .- 30*(x.^2) .+ 3)];
    fx = [v0 v1 3*x (7.5*(x.^2) .- 1.5) (1.0/8.0)*(140*(x.^3) .- 60*x)];
else
    error("polynomial degree must be less than 5");
end

fx = 2.0*fx;

return f, fx;

end
