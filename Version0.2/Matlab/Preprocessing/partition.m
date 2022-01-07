function [epart, npart] = partition(t2f,ne,np,metis)

% find metis executable
%metis = findexec(metis,version);

[metisstatus0,~] = system("which " + metis);
[metisstatus1,~] = system("which mpmetis");
[metisstatus2,~] = system("which /usr/bin/mpmetis");
[metisstatus3,~] = system("which /usr/local/bin/mpmetis");
[metisstatus4,~] = system("which /opt/local/bin/mpmetis");        

if metisstatus0==0
elseif metisstatus1==0        
    metis = "mpmetis";        
elseif metisstatus2==0        
    metis = "/usr/bin/mpmetis";    
elseif metisstatus3==0        
    metis = "/usr/local/bin/mpmetis";    
elseif metisstatus4==0        
    metis = "/opt/local/bin/mpmetis";    
else            
    error("Exasim search in /usr/bin, /usr/local/bin, and /opt/local/bin and could not find Metis. Please see the documentation to install it. After installation, please set its path to app.metis"); 
end

disp('Writing input files for METIS...');
dlmwrite('temp.txt', ne, 'delimiter', ' ','precision',10);
dlmwrite('temp.txt', t2f, '-append', 'delimiter', ' ','precision',10);

% call mpmetis
str = ['!' char(metis) ' temp.txt ' num2str(np)];
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










% function [epart, npart] = partition(t2f,ne,np,weightingElements,elcon)
% 
% % Expected formats:
% % t2f: ne x nve [2D array]
% % elcon: (npf*nfe) x ne [2D array]
% 
% if nargin<4; weightingElements = 0; end
% if nargin<5; elcon = []; weightingElements = 0; end
% 
% if size(t2f,2) == ne; t2f = t2f'; end
% if size(elcon,1) == ne; elcon = elcon'; end
% 
% % current directory
% cdir = pwd;
% 
% % move to directory that contains metis programs
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
%     cd(strcat(ps(1:is(end-up)),'metis/mac'));
% elseif isunix
%     cd(strcat(ps(1:is(end-up)),'metis/linux'));
% end
% 
% % Generate a temporary file to be used in METIS
% if weightingElements == 0 || isempty(elcon)
%     disp('Writing input files for METIS...');
%     dlmwrite('temp.txt', ne, 'delimiter', ' ','precision',10);
%     dlmwrite('temp.txt', t2f, '-append', 'delimiter', ' ','precision',10);
% else
%     disp('Computing element weights for METIS mesh partition...');
%     elcon_sorted = sort(elcon,1);
%     elcon_sorted = elcon_sorted';
%     elcon_sorted(:,end+1) = elcon_sorted(:,end);
%     elcon_sorted = elcon_sorted';
%     [a,b] = find((elcon_sorted(2:end,:)-elcon_sorted(1:end-1,:)) <= 0);
%     [~,i_tmp,~] = unique(b);
%     if length(i_tmp) ~= ne; error('Something wrong.'); end
%     weights = a(i_tmp).^2;     % We assume cost is O(ncf^2). Actual cost of matrix assembly is O(1) + O(ncf) + O(ncf^2). ncf is the number of unique traced nodes on the element faces
%     % Note: We take wieght = ncf^2 to balance the cost of matrix-vector product and
%     %       preconditioner solve among processors, as these operations are proportional to the square of
%     %       the number of nodes in the 0- and 1-level overlap node subdomain, respectively.
%     %       This is so since we make the domain decomposition for nodes such that (no. nodes
%     %       in 0-level overlap DD / no. elements in 0-level overlap DD) is approximtely constant for all processors
%     % Also, we note that if a few processors do much less work than the
%     %       average is much better than if a few processors do much more work
%     %       that the average.
%     
%     disp('Writing input files for METIS...');
%     dlmwrite('temp.txt', [ne, 1], 'delimiter', ' ','precision',10);
%     dlmwrite('temp.txt', [weights(:), t2f], '-append', 'delimiter', ' ','precision',10);
% end
% 
% disp('Calling METIS and reading output files...');
% if isempty(find(t2f(:,end)==-1, 1))==0  % fix -1
%     fin = fopen('temp.txt');
%     fout = fopen('output.txt','w');
%     while ~feof(fin)
%        s = fgetl(fin);
%        s = strrep(s, ' -1', '');              
%        fprintf(fout,'%s\n',s);       
%     end
%     fclose(fin);
%     fclose(fout);
%         
%     % call mpmetis
%     str = ['!./mpmetis output.txt ' num2str(np)];
%     eval(str);
% 
%     % get mesh partitioning data
%     str = ['output.txt.epart.' num2str(np)];
%     epart = textread(str,'%d');
% 
%     % get node partitioning data
%     str = ['output.txt.npart.' num2str(np)];
%     npart = textread(str,'%d');
% 
%     % remove files
%     delete('temp.txt');
%     delete('output.txt');
%     str = ['output.txt.epart.' num2str(np)];
%     delete(str);
%     str = ['output.txt.npart.' num2str(np)];
%     delete(str);
% 
%     % move back to current directory
%     cd(cdir);
% else
%     % call mpmetis
%     str = ['!./mpmetis temp.txt ' num2str(np)];
%     eval(str);
% 
%     % get mesh partitioning data
%     str = ['temp.txt.epart.' num2str(np)];
%     epart = textread(str,'%d');
% 
%     % get node partitioning data
%     str = ['temp.txt.npart.' num2str(np)];
%     npart = textread(str,'%d');
% 
%     % remove files
%     delete('temp.txt');
%     str = ['temp.txt.epart.' num2str(np)];
%     delete(str);
%     str = ['temp.txt.npart.' num2str(np)];
%     delete(str);
% 
%     % move back to current directory
%     cd(cdir);    
% end
