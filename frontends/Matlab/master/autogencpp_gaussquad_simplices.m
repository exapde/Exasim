
for nd = [2,3]
    
    if nd == 2
        pgauss_max = 16;
        fileName = 'gaussQuadTri.cpp';
        str1 = 'GAUSSQUAD_TRI';
        str3 = 'gaussQuadTri';
    elseif nd == 3
        pgauss_max = 15;
        fileName = 'gaussQuadTet.cpp';
        str1 = 'GAUSSQUAD_TET';
        str3 = 'gaussQuadTet';
    else
        error('Function only valid for nd == 2 or nd == 3');
    end

    delete(fileName);
    gid = fopen(fileName,'wt');

    str = ['#ifndef __',str1];
    fprintf(gid, '%s\n',str);
    str = ['#define __',str1];
    fprintf(gid, '%s\n',str);
    fprintf(gid, '\n');
    str = '// Written by: C. Nguyen & P. Fernandez';
    fprintf(gid, '%s\n',str);
    fprintf(gid, '\n');
    str = ['void ',str3,'(vector<double> *x_p, vector<double> *w_p, Int *nq, Int pgauss)'];
    fprintf(gid, '%s\n',str);
    str = '{';
    fprintf(gid, '%s\n',str);
    
    str = 'double *x;';
    fprintf(gid, '\t%s\n',str);
    str = 'double *w;';
    fprintf(gid, '\t%s\n',str);
    
    str = 'switch (pgauss) {';
    fprintf(gid, '\t%s\n',str);

    for pgauss=0:pgauss_max
        if nd == 2
            [x,w]=gaussquad2d(pgauss);
        elseif nd == 3
            [x,w]=gaussquad3d(pgauss);
        end

        str = ['case ',num2str(pgauss),':'];
        fprintf(gid, '\t\t%s\n',str);

        ng = length(w);
        str = ['*nq = ',num2str(ng),';'];
        fprintf(gid, '\t\t\t%s\n',str);
        str = ['x_p[0].resize(',num2str(ng*nd),');'];
        fprintf(gid, '\t\t\t%s\n',str);
        str = ['w_p[0].resize(',num2str(ng),');'];
        fprintf(gid, '\t\t\t%s\n',str);
        str = 'x = &x_p[0][0];';
        fprintf(gid, '\t\t\t%s\n',str);
        str = 'w = &w_p[0][0];';
        fprintf(gid, '\t\t\t%s\n',str);
        for i=1:ng*nd
            str = ['x[',num2str(i-1),'] = ',num2str(x(i),'%1.16e'),';'];
            fprintf(gid, '\t\t\t%s\n',str);
        end
        for i=1:ng
            str = ['w[',num2str(i-1),'] = ',num2str(w(i),'%1.16e'),';'];
            fprintf(gid, '\t\t\t%s\n',str);
        end
        str = 'break;';
        fprintf(gid, '\t\t\t%s\n',str);
    end

    str = 'default:';
    fprintf(gid, '\t\t%s\n',str);
    str = ['error("Only can handle pgauss <= ',num2str(pgauss_max),'.");'];
    fprintf(gid, '\t\t\t%s\n',str);

    str = '}';
    fprintf(gid, '\t%s\n',str);
    str = '}';
    fprintf(gid, '%s\n',str);
    fprintf(gid, '\n');
    str = '#endif';
    fprintf(gid, '%s\n',str);

    fclose(gid);
    
end
