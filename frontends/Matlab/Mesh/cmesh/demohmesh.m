Rc = 3;
th  = atan2(0.5,Rc);
cir = 2*pi/3;

[xf,yf] = read_foil('ht13foil');
xf = xf - 0.5;
yf = yf + Rc;

DR = sqrt(Rc^2+ 0.5^2) - Rc;

xfn = atan2(xf,yf);
yfn = sqrt(xf.^2 + yf.^2);
yfn = yfn - Rc - DR;
xfn = xfn/(2*th) + 0.5;
yfn = yfn/(2*th);

clf, plot(xfn,yfn,'-o'), axis equal, view(2);

[x,y] = hmeshparam;
figure(1);
surf(x, y, 0*x), axis equal, view(2);
[xm, ym] = hmeshmap(xfn, yfn, x, y, 6, 4);
clf, surf(xm, ym, 0*xn), axis equal, view(2); 

xmin = min(min(xm));
xmax = max(max(xm));
xm(:,1) = xmin;
xm(:,end) = xmax;
ymax = max(max(ym));
ym(end,:) = ymax;
ym(:,end) = ym(:,1);

xm = xm*2*th;
ym = ym*2*th;
ym = ym + Rc;
xm = xm-th;
xmin = min(min(xm));
xmax = max(max(xm));
ind = xm<-th;
xm(ind) = -th + (xm(ind)+th)*(-cir/2+th)/(xmin+th);
ind = xm>th;
xm(ind) =  th + (xm(ind)-th)*(cir/2-th)/(xmax-th);

xl = x;
yl = -y-0.001;
% surf(xl, yl, 0*xl), axis equal, view(2);
[xlm, ylm] = hmeshmap(xf, yf, xl, yl, 6, 4);
% pause;
% clf, surf(xlm, ylm, 0*xlm), axis equal, view(2);

%xm = (pi/6)*(xm-0.5)/0.75;
%ym = ym+4;

xn = ym.*sin(xm);
yn = ym.*cos(xm);

clf, surf(xn, yn, 0*xn), axis equal, view(2); 

% hold on;

% xm = [xlm];
% ym = [ylm];
% 
% %xm = (pi/6)*(xm-0.5)/0.75;
% %ym = ym+4;
% 
% xn = ym.*sin(xm);
% yn = ym.*cos(xm);
% 
% surf(xm, ym, 0*xm), axis equal, view(2);