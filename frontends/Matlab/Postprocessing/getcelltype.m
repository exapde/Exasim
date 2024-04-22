function [cell_t] = getcelltype(nd,elemtype)

if nd==2
    if elemtype==0
        cell_t = 5;
    else
        cell_t = 9;
    end
elseif nd==3
    if elemtype==0
        cell_t = 10;
    else
        cell_t = 12;
    end
else
    error('Number of mesh spatial dimensions should be 2 or 3')
end



