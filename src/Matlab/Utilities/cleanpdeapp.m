cdir = pwd;
ii = strfind(cdir, "Exasim"); 
testdir = cdir(1:(ii+5)) + "/tests";
cd(testdir);

mfiles = dir('**/*.m');
n = length(mfiles);
for i = 1:n      
    if string(mfiles(i).name) == "pdeapp.m"
        filein = string(mfiles(i).folder) + '/' + string(mfiles(i).name); 
        fid = fopen(filein);
        s = fscanf(fid,'%c');
        fclose(fid);        
        ii = strfind(s,'% visualize');
        if isempty(ii)==0
            ii = ii(1);
            s(ii:end) = [];                    
            fid = fopen(filein,'w');        
            fprintf(fid, '%c', s);
            fclose(fid);     
        end
        ii = strfind(s,'pde.dt'); 
        if isempty(ii)==0
            ii = ii(1);
            jj = strfind(s(ii:end),','); jj = jj(1);
            kk = strfind(s(ii:end),')'); kk = kk(1);
            s = strrep(s,s((ii+jj):(ii+kk-1)),'200)');        

            ii = strfind(s,'pde.soltime'); ii = ii(1);
            jj = strfind(s(ii:end),';'); jj = jj(1);
            s = strrep(s,s(ii:(ii+jj-1)),'pde.soltime = 10:10:200;');        

            fid = fopen(filein,'w');        
            fprintf(fid, '%c', s);
            fclose(fid);     
        end
    end    
end

jfiles = dir('**/*.jl');
n = length(jfiles);
for i = 1:n      
    if string(jfiles(i).name) == "pdeapp.jl"
        filein = string(jfiles(i).folder) + '/' + string(jfiles(i).name); 
        fid = fopen(filein);
        s = fscanf(fid,'%c');
        fclose(fid);        
        ii = strfind(s,'# visualize');
        if isempty(ii)==0
            ii = ii(1);
            s(ii:end) = [];                    
            fid = fopen(filein,'w');        
            fprintf(fid, '%c', s);
            fclose(fid);     
        end
        ii = strfind(s,'pde.dt'); 
        if isempty(ii)==0
            ii = ii(1);
            jj = strfind(s(ii:end),'('); jj = jj(1);
            kk = strfind(s(ii:end),')'); kk = kk(1);
            s = strrep(s,s((ii+jj):(ii+kk-1)),'200');        

            ii = strfind(s,'pde.soltime'); ii = ii(1);
            jj = strfind(s(ii:end),';'); jj = jj(1);
            s = strrep(s,s(ii:(ii+jj-1)),'pde.soltime = collect(10:10:length(pde.dt));');        

            fid = fopen(filein,'w');        
            fprintf(fid, '%c', s);
            fclose(fid);     
        end        
    end    
end


pfiles = dir('**/*.py');
n = length(pfiles);
for i = 1:n      
    if string(pfiles(i).name) == "pdeapp.py"
        filein = string(pfiles(i).folder) + '/' + string(pfiles(i).name); 
        fid = fopen(filein);
        i
        filein
        s = fscanf(fid,'%c');
        fclose(fid);        
        ii = strfind(s,'# visualize');
        if isempty(ii)==0
            ii = ii(1);
            s(ii:end) = [];                    
            fid = fopen(filein,'w');        
            fprintf(fid, '%c', s);
            fclose(fid);     
        end
        ii = strfind(s,"pde['dt']"); 
        if isempty(ii)==0
            ii = ii(1);
            jj = strfind(s(ii:end),'('); jj = jj(1);
            kk = strfind(s(ii:end),')'); kk = kk(1);
            s = strrep(s,s((ii+jj):(ii+kk-2)),'200');        

            ii = strfind(s,"pde['soltime']"); ii = ii(1);
            jj = strfind(s(ii:end),';'); jj = jj(1);
            s = strrep(s,s(ii:(ii+jj-1)),"pde['soltime'] = numpy.arange(10,pde['dt'].size+1,10);");        

            fid = fopen(filein,'w');        
            fprintf(fid, '%c', char(s));
            fclose(fid);     
        end                
    end    
end

cd(cdir);




