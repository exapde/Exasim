function face = getelemface(dim,elemtype)

switch dim
    case 1
        face=[1;2];
    case 2
        if elemtype==0
            face=[[2,3];[3,1];[1,2]];
        elseif elemtype==1
            face=[[1,2];[2,3];[3,4];[4,1]];
        end
    case 3
        if elemtype==0
            face=[[2,3,4];[1,4,3];[1,2,4];[1,3,2]];
        elseif elemtype==1
            face=[[1,4,3,2];[5,6,7,8];[1,2,6,5];[3,4,8,7];[2,3,7,6];[4,1,5,8]];
        end
    otherwise
        error('Only can handle dim=1, dim=2 or dim=3');
end
face = face';


