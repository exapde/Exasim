nt = size(sol,4);
h=figure;
shwmsh=0;
for n = 1:nt
    t = plotsol(n, 0);

    make_animation( h,n,'sample_animation.gif' )
    pause(0.05) %you can enter the time in pause to change the loop
end
function make_animation( h,index,filename )
drawnow
frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
if index == 1
    imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
else
    imwrite(imind,cm,filename,'gif','WriteMode','append');
end
end