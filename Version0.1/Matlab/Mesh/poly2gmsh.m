function poly2gmsh(fname, pv, elemtype)
    
    if nargin<3 
        elemtype = 0;
    end
        
    fid = fopen(fname, 'w');
    
    npv = size(pv,1);
    
    for i = 1:npv
        fprintf(fid, 'Point(%d) = {%g, %g, 0};\n', i, pv(i,:));
    end
    fprintf(fid, '\n');
    for i = 1:npv
        fprintf(fid, 'Line(%d) = {%d, %d};\n', i, i, mod(i,npv) + 1);
    end    
    fprintf(fid, '\n');

    str = int2str(1:npv);
    commapos = find(str(2:end)==' ' & str(1:end-1)~=' ') + 1;
    str(commapos) = ',';

    fprintf(fid, 'Line Loop(%d) = {%s};\n', npv, str);
    fprintf(fid, 'Plane Surface(1000) = {%d};\n', npv);    
    if elemtype==1
        fprintf(fid, 'Recombine Surface {1000};\n');
    end
    
    fclose(fid);

end



