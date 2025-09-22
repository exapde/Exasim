function [fnew,t2fnew,f2fnew,flev] = faceordering(f, t2f)

% number of faces
nf  = size(f,1);    

% get face-to-face connectivities
f2f = mkf2f(f, t2f);

% boundary faces are ordered first
fa   = find(f(:,end) < 0); 
fc   = fa;
flev = length(fc);

while 1
    % get face-to-face connectivities for faces fa
    f2fn = f2f(:,fa);        
    
    % get all neighboring faces of fa
    fb   = unique(f2fn(2:end,:)); 
    
    % remove zero
    if fb(1) == 0                 
        fb = fb(2:end);
    end 
    % obtain faces in fb that are not in [fa; fc]
    fim  = ismember(fb,[fa;fc]);
    im   = find(fim==0);
    
    % update fa
    fa   = fb(im);        
        
    % update fc and flev
    fc   = [fc; fa];
    flev = [flev length(fc)];
    
    if length(fc) == nf
        break;
    end
end

% obtain new face-to-element connectivities
fnew = f(fc,:);

% update new element-to-face connectivities
[nt nfv] = size(t2f);
a      = zeros(nf,1);
a(fc)  = (1:nf)';
t2f    = reshape(t2f,nfv*nt,1);
t2fnew = a(t2f);
t2fnew = reshape(t2fnew,nt,nfv);

% update f2f
f2fnew = mkf2f(fnew, t2fnew);


% nf  = size(f,1);    % number of faces
% 
% % boundary faces are ordered first
% fc   = find(f(:,end) < 0); 
% fd   = fc;
% flev = length(fd);
% 
% while 1
%     % get face-to-face connectivities for old faces fc
%     f2fn = f2f(:,fc);
%     
%     % current faces
%     fa   = f2fn(1,:)';      
%     
%     % get all neighboring faces of fa
%     fb   = unique(f2fn(2:end,:)); 
%     
%     % remove zero
%     if fb(1) == 0                 
%         fb = fb(2:end);
%     end 
%     % obtain faces in fb that are not in [fa; fd]
%     fim  = ismember(fb,[fa;fd]);
%     im   = find(fim==0);
%     fc   = fb(im);        
%         
%     % update fd and flev
%     fd   = [fd; fc];
%     flev = [flev length(fd)];
%     
%     if length(fd) == nf
%         break;
%     end
% end
% 
% fnew = f(fd,:);


%     % remove faces in fb that are in fd
%     fim  = ismember(fb,fd);
%     im   = find(fim==1);
%     fb(im) = [];        
%     
%     % remove faces in fb that are in fa  
%     fin  = ismember(fb,fa);  
%     in   = find(fin==0);
%     fc   = fb(in);            


