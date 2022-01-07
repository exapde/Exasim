function screen2png(filename,frame,sz,res)
% SCREEN2PNG Generate a PNG file of the current figure with
% dimensions consistent with the figure's screen dimensions.
%
% SCREEN2PNG('filename') saves the current figure to the
% PNG file "filename".
%

if nargin < 1
error('Not enough input arguments!')
end

if nargin<3, sz=[]; end
if nargin<4, res=100; end

fn=sprintf('%s%05d.png',filename,frame);

colormap(jet(65536));
delete(findobj(gcf,'tag','Colorbar'));
%axis equal,axis tight; 
axis off

oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');

if nargin<3 || isempty(sz)
    newpos = scrpos/100;
    newpos(end) = newpos(end)/1.25;
else
    newpos = sz;
end

set(gcf,'PaperUnits','inches','PaperPosition',newpos);
print('-dpng','-z',sprintf('-r%d',res),fn);

drawnow
set(gcf,'Units',oldscreenunits,...
'PaperUnits',oldpaperunits,...
'PaperPosition',oldpaperpos);
