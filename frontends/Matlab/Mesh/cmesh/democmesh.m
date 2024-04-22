[x,y] = cmeshparam4(45, 60, 60, 50, [40, 10, 10, 5, 1, 1, 5], [50, 1000, 200, 1000, 50]);
figure(1);
surf(x, y, 0*x), axis equal, view(2);

[xf,yf] = read_foil('ht13foil');

pause;
clf, plot(xf,yf,'-o'), axis equal, view(2);

[xm, ym] = cmeshmap(xf, yf, x, y, 6, 4);

pause;
clf, surf(xm, ym, 0*xm), axis equal, view(2);

[xf,yf] = read_foil('sd7003foil');

pause;
clf, plot(xf,yf,'-o'), axis equal, view(2);

[xm, ym] = cmeshmap(xf, yf, x, y, 6, 4);

pause;
clf, surf(xm, ym, 0*xm), axis equal, view(2);


porder=4;
[x,y] = cmeshparam6(11*porder+1, 7*porder+1, 7*porder+1, 7*porder+1, 7*porder+1, 12*porder+1, ...
                    [100, 20, 20, 20, 20, 1, 1, 1, 1, 1, 1], [5, 100, 50, 20, 50, 100, 5]*10);
                
[xm, ym] = cmeshmap(xf, yf, x, y, 6, 4);

pause;
clf, surf(xm, ym, 0*xm), axis equal, view(2);
