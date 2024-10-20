function poly2gmsh(fname, pv, h)
    
    fid = fopen(fname, 'w');
    
    npv = size(pv,1);
    
    if nargin>2
        for i = 1:npv
            fprintf(fid, 'Point(%d) = {%g, %g, 0, %g};\n', i, [pv(i,:) h]);
        end        
    else
        for i = 1:npv
            fprintf(fid, 'Point(%d) = {%g, %g, 0};\n', i, pv(i,:));
        end
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
    fprintf(fid, 'Plane Surface(1) = {%d};\n', npv);
    fprintf(fid, 'Recombine Surface {1};\n');    

    fclose(fid);

end
