function plotboundary(mesh)

nb = max(mesh.f(:));
colors = lines(nb);
figure(1); clf; boundaryplot(mesh,1,colors(1,:));
hold on;
for i = 2:nb
  boundaryplot(mesh,i,colors(i,:));
end
