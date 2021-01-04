function xiny(x,y,opt::Int=1)
# Determine if each row of x is a member of y
# If row j of x is a member of y and x(j,:) = y(k,:) then in[j] = k
# Else in[j] = 0

if opt == 1
    m,dim = size(x);
else
    dim,m = size(x);
end

in = zeros(Int,m,1);

if opt==1
    if dim==1
        for j=1:m
            d2 = (y[:,1].-x[j,1]).^2;
            md,id = findmin(d2);
            if md<1e-12
                in[j]=id;
            end
        end
    elseif dim==2
        for j=1:m
            d2 = (y[:,1].-x[j,1]).^2 + (y[:,2].-x[j,2]).^2;
            md,id = findmin(d2);
            if md<1e-12
                in[j]=id;
            end
        end
    elseif dim==3
        for j=1:m
            d2 = (y[:,1].-x[j,1]).^2 + (y[:,2].-x[j,2]).^2 + (y[:,3].-x[j,3]).^2;
            md,id = findmin(d2);
            if md<1e-12
                in[j]=id;
            end
        end
    else
        n = size(y,1);
        for j=1:m
            d2 = sum((y - repmat(x[j,:],[n 1])).^2,2);
            md,id = findmin(d2);
            if md<1e-12
                in[j]=id;
            end
        end
    end
else
    if dim==1
        for j=1:m
            d2 = (y[1,:].-x[1,j]).^2;
            md,id = findmin(d2);
            if md<1e-12
                in[j]=id;
            end
        end
    elseif dim==2
        for j=1:m
            d2 = (y[1,:].-x[1,j]).^2 + (y[2,:].-x[2,j]).^2;
            md,id = findmin(d2);
            if md<1e-12
                in[j]=id;
            end
        end
    elseif dim==3
        for j=1:m
            d2 = (y[1,:].-x[1,j]).^2 + (y[2,:].-x[2,j]).^2 + (y[3,:].-x[3,j]).^2;
            md,id = findmin(d2);
            if md<1e-12
                in[j]=id;
            end
        end
    else
        n = size(y,2);
        for j=1:m
            d2 = sum((y - repmat(x[:,j],[1 n])).^2,1);
            md,id = findmin(d2);
            if md<1e-12
                in[j]=id;
            end
        end
    end
end

return in;

end
