function getelemface(dim,elemtype)

if dim== 1
    face=[1;2];
    nfe = 2;
    nvf = 1;
elseif dim==2
    if elemtype==0
        face=[[2,3];[3,1];[1,2]];
        nfe = 3;
        nvf = 2;
    elseif elemtype==1
        face=[[1,2];[2,3];[3,4];[4,1]];
        nfe = 4;
        nvf = 2;
    end
elseif dim==3
    if elemtype==0
        face=[[2,3,4];[1,4,3];[1,2,4];[1,3,2]];
        nfe = 4;
        nvf = 3;
    elseif elemtype==1
        face=[[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]];
        nfe = 6;
        nvf = 4;
    end
else
    error("Only can handle dim=1, dim=2 or dim=3");
end

face = reshape(face,nvf,nfe);

return face;

end
