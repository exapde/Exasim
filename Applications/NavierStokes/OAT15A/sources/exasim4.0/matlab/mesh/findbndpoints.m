function ip = findbndpoints(mesh,bnd)
% ip = findbndpoints(mesh,bnd)
% Function to return the index numbers of all points on boundary bnd
% in a given HDG mesh. Points are returned in a geometrically sequential
% order, starting at the first point on the first face of the mesh
% which is on boundary bnd.
%
% First and last points in output "ip" are repeated.
%
%--------------------------------------------------------------------------
% REVISION HISTORY
% When     Who               What
% 08Nov12  Hemant Chaurasia  Created
%--------------------------------------------------------------------------

bf = find(mesh.f(:,4)==-bnd);
nbf = size(bf,1);

ip = mesh.f(bf(1),1:2);        % 2 points from first face on bnd, natural order
iface = bf(1);
fdone = 1;
while fdone<nbf
    nextf = bf(find(mesh.f(bf,1)==mesh.f(iface(end),2)));
    if length(nextf)==1
        iface(end+1) = nextf;
    else
        error('Cannot be applied to boundaries around zero-thickness objects.');
    end
    ip(end+1) = mesh.f(iface(end),2);
    fdone = fdone + 1;
end