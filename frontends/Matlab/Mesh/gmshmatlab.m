function gmshmatlab(filename,str)

% % current directory
% cdir = pwd;

% % move to directory that contains gmsh program
% if ispc
%     sslash = '\';
% elseif isunix
%     sslash = '/';
% end
% ps=strcat(pwd,sslash);
% is=find(ps==sslash);
% up=0;
% while ~(strcmp(ps(is(end-up-1)+1:is(end-up)-1),'Exasim') || strcmp(ps(is(end-up-1)+1:is(end-up)-1),'EXASIM'))
%     up = up+1;
% end
% if ismac
%     mdir = strcat(ps(1:is(end-up)),'gmsh');
% elseif isunix
%     mdir = strcat(ps(1:is(end-up)),'gmsh');
% end
% cd(mdir);

% % copy geo files to gmsh directory
% cfile = strcat(cdir,[sslash '*.geo']);
% copyfile(cfile,mdir);

% call gmsh
str = ['!gmsh ' filename '.geo ' str];
eval(str);

% % copy geo file to gmsh directory
% cfile = strcat(mdir,[sslash filename '.msh']);
% copyfile(cfile,cdir);
% 
% % move back to current directory
% cd(cdir);    



