function [elem2cpu,extintelem,extintelempts,face2cpu,extintface,extintfacepts,ent2cpu,extintent,extintentpts,outelem] =...
         overlappingpartition(t2f,t,elcon,f,hybrid,overlappinglevel,nproc,check,ratioMaxToAvgEntPerProcessor,weightingDD) 
     
if nargin < 8; check = 0; end
if nargin < 9; ratioMaxToAvgEntPerProcessor = 1.0e6; end
if nargin < 10; weightingDD = 0; end

ne = size(t,2);
% nf = size(facecon,2);
% npf = size(facecon,1);
% nfe = size(t2f,1);

% Nonoverlapping element and face partition
disp(' '); disp('Calling t2fpartition...');
[elem2cpu, face2cpu] = t2fpartition(t2f+1,ne,nproc,weightingDD,elcon);    

% % % % % % % % % % % % % Report quality of METIS domain decomposition for faces and elements:
% % % % % % % % % % % % elem2cpu_sorted = sort(elem2cpu);
% % % % % % % % % % % % elem2cpu_sorted(end+1) = elem2cpu_sorted(end);
% % % % % % % % % % % % [a,b] = find((elem2cpu_sorted(2:end)-elem2cpu_sorted(1:end-1)) == 0);
% % % % % % % % % % % % [~,i_tmp,~] = unique(b);
% % % % % % % % % % % % if length(i_tmp) ~= nproc; error('Something wrong.'); end
% % % % % % % % % % % % numElemPerProc = a(i_tmp);
% % % % % % % % % % % % disp(['Minimum elements per processor: ', num2str(min(numElemPerProc))]);
% % % % % % % % % % % % disp(['Maximum elements per processor: ', num2str(max(numElemPerProc))]);
% % % % % % % % % % % % disp(['Average elements per processor: ', num2str(ne / nproc)]);

% Overlapping HDG patition
disp(' '); disp('Calling hdgpartition...');
[extintelem,extintelempts,extintface,extintfacepts,outelem] = hdgpartition(t,t2f,elem2cpu,face2cpu,overlappinglevel,nproc);
% [extintelemc,extintelemptsc,extintfacec,extintfaceptsc,outelemc] = hdgpartitionc(t,t2f,elem2cpu,face2cpu,[size(t) nf nfe overlappinglevel nproc]);                         
% for i = 1:nproc
%     [max(abs(extintelem{i}(:)-extintelemc{i}(:))) max(abs(extintface{i}(:)-extintfacec{i}(:))) max(abs(outelemc{i}(:)-outelemc{i}(:)))]
%     [max(abs(extintelempts(:)-extintelemptsc(:))) max(abs(extintfacepts(:)-extintfaceptsc(:)))]
% end

% Overlapping EDG patition
if strcmp(hybrid,'hdg')==0
    %ndh = max(facecon(:))+1;
    %[ent2cpu,extintent,extintentpts] = edgpartitionc(extintelem, elcon, facecon, t2f, face2cpu, [ne nf ndh nfe npf nproc]);    
    disp(' '); disp('Calling mkv2t...'); disp(' ');
    [re,ce] = mkv2t(t+1,ne);
    disp(' '); disp('Calling edgpartition...');
    ent2cpu = edgpartition(elcon, f, t, t2f, re, ce, elem2cpu, face2cpu, nproc, overlappinglevel, ratioMaxToAvgEntPerProcessor);            
    extintent = cell(nproc,1);
    extintentpts = zeros(nproc,2);
    disp(' '); disp('Computing extintent and extintentpts...');
    for i = 1:nproc
        intent = find(ent2cpu==(i-1));            
        ent = elcon(:,extintelem{i}(1:sum(extintelempts(i,1:2)))+1)+1;
        ent = unique(ent(:));                                
        ent = ent(ent>0);                      
        extent = setdiff(ent,intent);      
        extintent{i}=[intent; extent]-1;
        extintentpts(i,:)=[length(intent) length(extent)];                       
    end    
else
    ent2cpu = []; extintent = []; extintentpts = [];
end

% TODO: Reorder elements and entities

for i = 1:nproc    
    if isempty(outelem{i})
        error('HDG element partition is incorrect.');
    end
    if ~isempty(intersect(extintelem{i},outelem{i}))
        error('HDG element partition is incorrect.');
    end
    extintelem{i} = [extintelem{i}; outelem{i}];
    extintelempts(i,3) = length(outelem{i});
end

% Check METIS decomposition
disp(' '); disp('Checking METIS decomposition...');
checkElements = 0;
checkEntities = 0;
for i = 1:nproc                           
    intelem = find(elem2cpu==(i-1));     
    checkElements = checkElements + length(intelem);
    
    intface = find(face2cpu==(i-1));     
    checkEntities = checkEntities + length(intface);
    
    face = t2f(:,intelem)+1;  
    face = unique(face(:));
    if ~isempty(setdiff(intface,face)) 
        error('Some faces in the processor are not connected to any element in the processor.');
    end                  
    
    n = size(f,1);
    tm = f((n-1):n,intface)+1;  
    tm = unique(tm(:));
    tm = tm(tm>0);
    if ~isempty(setdiff(intelem,tm)) 
        error('Some elements in the processor are not connected to any face in the processor.');
    end                            
end

numElem = ne;
numEntities = max(t2f(:))+1;
if numEntities ~= length(unique(t2f(:))); error('It looks like there are ghost faces in the mesh.'); end
if checkElements ~= numElem; error('Domain decomposition for elements was not performed properly.'); end
if checkEntities ~= numEntities; error('Domain decomposition for entities was not performed properly.'); end

% % Check that all elements in the processor have at least one entity in the
% % processor
% disp(' ');
% disp('Check that all elements in the processor have at least one entity in the processor...');
% disp(' ');
% for i=1:nproc
%     disp(['Proc. No. ', num2str(i),'/',num2str(nproc)]);
%     elementsToRemove = [];
%     for j=1:length(intelem{i})
%         if length(setdiff(intface{i},t2f(:,intelem{i}(j))+1)) == length(intface{i})
%             entitiesInElement = t2f(:,intelem{i}(j))+1;
%             candidateProcessors = entpart(entitiesInElement);
%             k = mode(candidateProcessors);
%             intelem{k} = sort([intelem{k};intelem{i}(j)]);
%             elementsToRemove = [elementsToRemove;intelem{i}(j)];
%             error('An element in the processor had did not have any entity in the processors. The issue has been fixed successfully.');
%         end
%     end
% %     if ~isempty(elementsToRemove)
% %         intelem{i} = sort(setdiff(intelem{i},elementsToRemove));
% %     end
% end

if check==1
    % check HDG patition
    [re,ce] = mkv2t(t+1,ne);                         
    for i = 1:nproc
        intelem = find(elem2cpu==(i-1));     

        intface = find(face2cpu==(i-1));     
        
        elem = intelem;
        for j=1:overlappinglevel            
            elem = node2elem(t(:,elem)+1,re,ce);
        end
        extelem = setdiff(elem,intelem);                        
        if max(abs(extintelem{i}(1:sum(extintelempts(i,1:2)))+1-[intelem; extelem]))~=0
            error('HDG element partition is incorrect.');
        end
        if max(abs(extintelempts(i,1:2)-[length(intelem) length(extelem)]))~=0
            error('HDG element partition is incorrect.');
        end
        elem = intelem;
        for j=1:overlappinglevel+1            
            elem = node2elem(t(:,elem)+1,re,ce);
        end
        extelem = setdiff(elem,intelem);     
        if max(abs(outelem{i}+1-setdiff([intelem; extelem],extintelem{i}(1:sum(extintelempts(i,1:2)))+1)))~=0            
            error('HDG element partition is incorrect.');
        end
        
        face = t2f(:,extintelem{i}(1:sum(extintelempts(i,1:2)))+1)+1;  
        face = unique(face(:));
        extface = setdiff(face,intface);                        
        if max(abs(extintface{i}+1-[intface; extface]))~=0
            error('HDG face partition is incorrect.');
        end
        if max(abs(extintfacepts(i,:)-[length(intface) length(extface)]))~=0
            error('HDG face partition is incorrect.');
        end
    end

%     if strcmp(hybrid,'hdg')==0                                 
%         % check EDG patition
%         ent2cpum = edgpartition(facecon, t, t2f, re, ce, elem2cpu, face2cpu, nproc, overlappinglevel);        
%         if max(abs(ent2cpu-ent2cpum))~=0
%             error('EDG entity partition is incorrect.');
%         end
%         for i = 1:nproc
%             intent = find(ent2cpum==(i-1));            
%             ent = elcon(:,extintelem{i}(1:sum(extintelempts(i,1:2)))+1)+1;
%             ent = unique(ent(:));                                
%             ent = ent(ent>0);                      
%             extent = setdiff(ent,intent);      
%             if max(abs(extintent{i}(1:extintentpts(i,1))+1-intent))~=0
%                 error('EDG entity partition is incorrect.');
%             end
%             if max(abs(extintent{i}+1-[intent; extent]))~=0
%                 error('EDG entity partition is incorrect.');
%             end
%             if max(abs(extintentpts(i,:)-[length(intent) length(extent)]))~=0
%                 error('EDG entity partition is incorrect.');
%             end            
%         end
%     end
end

% elempart = elempart+1;
% entpart = entpart+1;
% 
% % loop over each subdomain
% checkElements = 0;
% checkEntities = 0;
% for i=1:nproc    
%     % list of elements in subdomain i                       
%     intelem{i} = find(elempart==i); 
%     checkElements = checkElements + length(intelem{i});
%     
%     % list of edg nodes in subdomain i                       
%     indent{i} = find(entpart==i); 
%     checkEntities = checkEntities + length(indent{i});
% end
% 
% numElem = mesh.ne;
% numEntities = max(mesh.t2f(:));
% if numEntities ~= length(unique(mesh.t2f(:))); error('It looks like there are ghost entities in the mesh.'); end
% if checkElements ~= numElem; error('Domain decomposition for elements was not performed properly.'); end
% if checkEntities ~= numEntities; error('Domain decomposition for entities was not performed properly.'); end
% 
% % Check that all elements in the processor have at least one entity in the
% % processor
% disp(' ');
% disp('Check that all elements in the processor have at least one entity in the processor...');
% disp(' ');
% for i=1:nproc
%     disp(['Proc. No. ', num2str(i),'/',num2str(nproc)]);
%     elementsToRemove = [];
%     for j=1:length(intelem{i})
%         if length(setdiff(indent{i},mesh.t2f(intelem{i}(j),:))) == length(indent{i})
%             entitiesInElement = mesh.t2f(intelem{i}(j),:);
%             candidateProcessors = entpart(entitiesInElement);
%             k = mode(candidateProcessors);
%             intelem{k} = sort([intelem{k};intelem{i}(j)]);
%             elementsToRemove = [elementsToRemove;intelem{i}(j)];
%             warning('An element in the processor had did not have any entity in the processors. The issue has been fixed successfully.');
%         end
%     end
%     if ~isempty(elementsToRemove)
%         intelem{i} = sort(setdiff(intelem{i},elementsToRemove));
%     end
% end
