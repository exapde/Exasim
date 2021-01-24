function UDG = initu(mesh,param)
%INITQ Initialize vector of unknowns
%    Q=INITU(MESH,APP,PARAM)
%
%    MESH:                    Mesh structure
%    APP:                     Application structure
%    PARAM{APP.NC;APP.NC}:    Cell array containing
%                             When param{i} is a double U(:,APP.NC,j,:) = param{j,i}
%                             When param{i} is a pointer to a function,
%                                           A(:,APP.NC,j,:) = param{j,i}(mesh.dgnodes)
%    Q(NPL,APP.NC,2,NT):      Scalar fucntion to be plotted
%

param = reshape(param',numel(param),[]);
UDG = zeros(size(mesh.dgnodes,1),length(param),size(mesh.dgnodes,3));
for i=1:length(param)
    if isa(param{i},'double')
        UDG(:,i,:) = param{i};
    else
        UDG(:,i,:) = param{i}(mesh.dgnodes);
    end
end
