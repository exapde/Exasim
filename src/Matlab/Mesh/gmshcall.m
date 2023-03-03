function [p,t] = gmshcall(pde, filename, nd, elemtype)

    % Mod to account for previously generated mesh - below lines are commented out
    
    % opts = "-format msh3";
    
    % % find gmsh executable
    % gmsh = findexec(pde.gmsh, pde.version);
    
    % disp('Gmsh mesh generator...');
    
    % % call gmsh
    % mystr = ['!' char(gmsh) ' ' char(filename) '.geo -' char(num2str(nd)) ' ' char(opts)];
    % eval(mystr);
    
    % [p,t] = gmsh2pt([char(filename) '.msh'],nd,elemtype);
    
    [p,t] = gmsh2pt([char(filename)],nd,elemtype);
    
    end
    
    function [p,t] = gmsh2pt(fname,nd,elemtype)
    
    fid = fopen(fname,'r');
    if fid == -1
        error('Cannot open file');
    end
    
    % 1	2-node line	Edge Lagrange P1
    % 2	3-node triangle	Triangle Lagrange P1
    % 3	4-node quadrangle	Quadrilateral Lagrange P1
    % 4	4-node tetrahedron	Tetrahedron Lagrange P1
    % 5	8-node hexahedron	Hexahedron Lagrange P1
    % 6	6-node prism	Wedge Lagrange P1    
    if nd==1 % line
        nve = 2;
        wcase = 1;
    elseif nd==2 && elemtype==0 % tri
        nve = 3;
        wcase = 2;
    elseif nd==2 && elemtype==1 % quad
        nve = 4;
        wcase = 3;
    elseif nd==3 && elemtype==0 % tet
        nve = 4;
        wcase = 4;
    elseif nd==3 && elemtype==1 % hex
        nve = 8;
        wcase = 5;
    end
    
    readuntil(fid, '$Nodes');
    np = fscanf(fid, '%d', 1);
    p = zeros(nd,np);
    for ii = 1:np
        foo = fscanf(fid, '%d', 1);
        p(1:nd,ii) = fscanf(fid, '%f', nd);
        fgetl(fid);
    end
    
    readuntil(fid, '$Elements');
    ne = fscanf(fid, '%d', 1);
    t = zeros(nve,ne);
    m = 0;
    for ii = 1:ne  
        eltype = fscanf(fid, '%d', 2);
        if (eltype(2) == wcase)
            m = m + 1;
            tii = fscanf(fid, '%d', 2+nve);
            t(:,m) = tii(3:end);
        end
        fgetl(fid);
    end
    t = t(:,1:m);
    fclose(fid);
    
    % nve = 0;
    % wcase = 0;
    % for ii = 1:ne    
    %     eltype = fscanf(fid, '%d', 2);
    %     if eltype(2) == 2
    %         nve = 3;
    %         tmp = fscanf(fid, '%d', 2+nve);
    %         t1 = tmp(3:end);        
    %     elseif (eltype(2) == 3) || (eltype(2) == 4)
    %         nve = 4;
    %         tmp = fscanf(fid, '%d', 2+nve);
    %         t1 = tmp(3:end);        
    %     elseif eltype(2) == 5
    %         nve = 8;
    %         tmp = fscanf(fid, '%d', 2+nve);
    %         t1 = tmp(3:end);                        
    %     end
    %     fgetl(fid);
    %     if nve>0
    %         wcase = eltype(2);
    %         break;
    %     end
    % end
    % 
    % if nve>0
    %     t = zeros(nve,ne);
    %     t(:,1) = t1;
    %     j = ii+1;
    %     for ii = j:ne
    %         eltype = fscanf(fid, '%d', 2);
    %         if eltype(2)==wcase
    %             tmp = fscanf(fid, '%d', 2+nve);
    %             t(:,ii) = tmp(3:end);      
    %         end
    %         fgetl(fid);
    %     end    
    %     t = t(:,any(t~=0,1));
    % else
    %     fclose(fid);
    %     error("Reading gmsh file failed");
    % end
    % fclose(fid);
    
    end
    
    function readuntil(fid,str)
        
        while ~feof(fid)
            fline = fgetl(fid);
            if ~isempty(fline) && strcmp(str, fline)
                break;
            end
        end
    end