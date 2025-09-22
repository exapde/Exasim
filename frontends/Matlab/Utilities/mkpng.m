function mkpng(mesh, sol, n, str, Minf)

% Create an invisible figure
f = figure('Visible','off');
ax = axes('Parent', f);

% Plot into the invisible axes
scaplot(mesh, eulereval(sol(:,:,:,n), str, 1.4, Minf), [], 1);
colormap(ax, jet);
set(ax, 'FontSize', 16);
axis(ax, 'off');
box(ax, 'on');
axis(ax, 'tight');

% Generate filename
fn = str + sprintf('%05d', n) + ".png";

% Export to file
exportgraphics(ax, fn, 'Resolution', 200);

% Clean up
close(f);

% function mkpng(fn,frame,res)
% 
% if nargin<3, res=100; end
% 
% axis equal,axis tight
% set(gcf,'Units','pixels');
% scrpos = get(gcf,'Position')
% sz = scrpos/100;
% 
% colormap(jet(65536));
% delete(findobj(gcf,'tag','Colorbar'));
% axis off
% fn=sprintf('%s%05d.png',fn,frame);
% set(gcf,'PaperUnits','inches','paperpos',sz,'papersize',[sz(3) sz(4)]);
% print('-dpng','-z',sprintf('-r%d',res),fn);
% 
% drawnow
% set(gcf,'Units',oldscreenunits,...
% 'PaperUnits',oldpaperunits,...
% 'PaperPosition',oldpaperpos)
