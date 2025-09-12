function [epart, npart] = t2fpartition(t2f,ne,np)
% t2f: ne * nve

sz = size(t2f);
if sz(2)==ne
    t2f = t2f';
end

% current directory
cdir = pwd;

% move to directory that contains metis programs
if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end
ps=strcat(pwd,sslash);
is=find(ps==sslash);
up=0;
while ~(strcmp(ps(is(end-up-1)+1:is(end-up)-1),'hdgv1.0') || strcmp(ps(is(end-up-1)+1:is(end-up)-1),'HDGv1.0'))
    up = up+1;
end
cd(strcat(ps(1:is(end-up)),'metis'));

% generate a temporary file to be used in metis
dlmwrite('temp.txt', ne, 'delimiter', ' ','precision',10);
dlmwrite('temp.txt', t2f, '-append', 'delimiter', ' ','precision',10);

% fix -1
if isempty(find(t2f(:,end)==-1, 1))==0
    fin = fopen('temp.txt');
    fout = fopen('output.txt','w');
    while ~feof(fin)
       s = fgetl(fin);
       s = strrep(s, ' -1', '');              
       fprintf(fout,'%s\n',s);       
    end
    fclose(fin);
    fclose(fout);
        
    % call mpmetis
    str = ['!./mpmetis output.txt ' num2str(np)];
    eval(str);

    % get mesh partitioning data
    str = ['output.txt.epart.' num2str(np)];
    epart = textread(str,'%d');

    % get node partitioning data
    str = ['output.txt.npart.' num2str(np)];
    npart = textread(str,'%d');

    % remove files
    delete('temp.txt');
    delete('output.txt');
    str = ['output.txt.epart.' num2str(np)];
    delete(str);
    str = ['output.txt.npart.' num2str(np)];
    delete(str);

    % move back to current directory
    cd(cdir);
else
    % call mpmetis
    str = ['!./mpmetis temp.txt ' num2str(np)];
    eval(str);

    % get mesh partitioning data
    str = ['temp.txt.epart.' num2str(np)];
    epart = textread(str,'%d');

    % get node partitioning data
    str = ['temp.txt.npart.' num2str(np)];
    npart = textread(str,'%d');

    % remove files
    delete('temp.txt');
    str = ['temp.txt.epart.' num2str(np)];
    delete(str);
    str = ['temp.txt.npart.' num2str(np)];
    delete(str);

    % move back to current directory
    cd(cdir);    
end


