function mkpng(fn,frame,res)

if nargin<3, res=100; end

axis equal,axis tight
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position')
sz = scrpos/100;

colormap(jet(65536));
delete(findobj(gcf,'tag','Colorbar'));
axis off
fn=sprintf('%s%05d.png',fn,frame);
set(gcf,'PaperUnits','inches','paperpos',sz,'papersize',[sz(3) sz(4)]);
print('-dpng','-z',sprintf('-r%d',res),fn);

drawnow
set(gcf,'Units',oldscreenunits,...
'PaperUnits',oldpaperunits,...
'PaperPosition',oldpaperpos)
