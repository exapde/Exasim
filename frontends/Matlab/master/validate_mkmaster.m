
clc

nodetype = 0;

for nd=1:3
    for porder=0:6
        for elemtype=0:1
            for quadOrder=1:3
                if quadOrder*porder <= 16
                    pgauss = quadOrder*porder;
                    
                    mesh.nd = nd;
                    mesh.porder = porder;
                    mesh.perm = [];
                    mesh.elemtype = elemtype;
                    [mesh.plocal,~,mesh.plocfc,~,~,~,~] = masternodes(porder,nd,elemtype,nodetype);
                    
                    master = mkmaster(mesh,pgauss);
                    
                    filename = ['P', num2str(porder), 'D', num2str(nd), 'E', num2str(elemtype), 'Q', num2str(quadOrder), '.bin'];

                    fileID = fopen(filename,'r');
                    data = fread(fileID,'double');    
                    fclose(fileID);
                    
                    start = 1;
                    N = length(master.plocvl(:));
                    endd = start + N - 1;
                    plocvl = data(start:endd);
                    errorNorm = norm(plocvl(:) - master.plocvl(:));
                    if errorNorm > 1.0e-9; disp([filename,' plocvl']); end
                        
                    start = endd + 1;
                    N = length(master.plocfc(:));
                    endd = start + N - 1;
                    plocfc = data(start:endd);
                    errorNorm = norm(plocfc(:) - master.plocfc(:));
                    if errorNorm > 1.0e-9; disp([filename,' plocfc']); end
                    
                    start = endd + 1;
                    N = length(master.gpvl(:));
                    endd = start + N - 1;
                    gpvlR = data(start:endd);
                    errorNorm = norm(gpvlR(:) - master.gpvl(:));
                    if errorNorm > 1.0e-9; disp([filename,' gpvl']); end
                    
                    start = endd + 1;
                    N = length(master.gwvl(:));
                    endd = start + N - 1;
                    gwvlR = data(start:endd);
                    errorNorm = norm(gwvlR(:) - master.gwvl(:));
                    if errorNorm > 1.0e-9; disp([filename,' gwvl']); end
                    
                    start = endd + 1;
                    N = length(master.gpfc(:));
                    endd = start + N - 1;
                    gpfcR = data(start:endd);
                    errorNorm = norm(gpfcR(:) - master.gpfc(:));
                    if errorNorm > 1.0e-9; disp([filename,' gpfc']); end
                    
                    start = endd + 1;
                    N = length(master.gwfc(:));
                    endd = start + N - 1;
                    gwfcR = data(start:endd);
                    errorNorm = norm(gwfcR(:) - master.gwfc(:));
                    if errorNorm > 1.0e-9; disp([filename,' gwfc']); end
                    
                    start = endd + 1;
                    N = length(master.shapvl(:));
                    endd = start + N - 1;
                    shapvlR = data(start:endd);
                    errorNorm = norm(shapvlR(:) - master.shapvl(:));
                    if errorNorm > 1.0e-9
                        disp([filename,' shapvl']);
                    end
                    
                    start = endd + 1;
                    N = length(master.shapfc(:));
                    endd = start + N - 1;
                    shapfcR = data(start:endd);
                    errorNorm = norm(shapfcR(:) - master.shapfc(:));
                    if errorNorm > 1.0e-9; disp([filename,' shapfc']); end
                end
            end
        end
    end
end
    