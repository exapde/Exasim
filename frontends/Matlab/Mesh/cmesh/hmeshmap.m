function [xgf,ygf] = hmeshmap( xf, yf, xr, yr, lw, ll)

small = 1.e-8;
shift = 0.04;                   % Good for SD7003

%nx = 41;
%ny = 21;

%[xr,yr] = read_foil();

dxu = xf(2)-xf(1);
dyu = yf(2)-yf(1);
dxl = xf(end-1)-xf(end);
dyl = yf(end-1)-yf(end);
tau = acos((dxu*dxl+dyu*dyl)/sqrt(dxu*dxu+dyu*dyu)/sqrt(dxl*dxl+dyl*dyl));
n   = 2-tau/pi;

chord = max(xf)-min(xf);
xmo = 0.5*(max(xf)+min(xf));
ymo = yf(1);
xn = 2.0*(1.0+shift)*(xf-xmo)/chord-shift;
yn = 2.0*(1.0+shift)*(yf-ymo)/chord;

% Near circle
[x,y] = trefftz_inv( xn, yn, n, 1/n, 0,1);

% plot(x,y,'-o'), axis equal, view(2);

% Center
xm = 0.5*(max(x)+min(x));
ym = 0.5*(max(y)+min(y));
fc = (max(x)-min(x)+max(y)-min(y))/4;
x = (x-xm)/fc;
y = (y-ym)/fc;

% Now a circle
[A,B] = GetTG( 200, x, y);

% Start mesh generation

% ix = xr(1,:)>1;
% xr(:,ix) = (xr(:,ix)-1)*(sqrt(lw+1)-1)+1;
% ix = xr(1,:)<-1;
% xr(:,ix) = (xr(:,ix)+1)*(sqrt(lw+1)-1)-1;
% yr = yr*sqrt(ll)+small;
% ix = (xr < 0);
% xr(ix) = -xr(ix);
% yr(ix) = -yr(ix);


%  Transformation
zr = complex(xr,yr);
wr = zr;

% Now a circle
[xc,yc] = trefftz_inv( 2.0*(real(wr)-0.5), 2.0*imag(wr), 2.0, 0.5, 0, 0);

surf(xc, yc, 0*xc);
view(2), axis equal;

% Now a near circle
[xg, yg] = TG( 2.0*xc, 2.0*yc, A, B);

% Finally the real thing
[xgf, ygf] = trefftz( xg*fc+xm, yg*fc+ym, n, 1/n, 0);

xgf = (xgf+shift)*chord/(2.0*(1.0+shift)) + xmo;
ygf = ygf*chord/(2.0*(1.0+shift)) + ymo;




function [x,y] = trefftz( x1, y1, n, cx, cy)
z1 = complex( x1, y1);
cc = complex( cx, cy);
A = ((z1-cc)./(z1+cc)).^n;
z = ((1+A)./(1-A))*n*cc;
x = real(z);
y = imag(z);

function [x,y] = trefftz_inv( x1, y1, n, cx, cy, track)
z1 = complex( x1, y1);
cc = complex( cx, cy);
A = ((z1-n*cc)./(z1+n*cc));
if track
   R = abs(A); 
   T = angle(A);
   for k = 2:size(T,1)
       d(1) = T(k) + 2*pi;
       d(2) = T(k);
       d(3) = T(k) - 2*pi;
       [minv,j] = min(abs(d-T(k-1)));
       T(k) = d(j);
   end
   B = ((R).^(1/n)).*exp(i*T/n);
else
   B = A.^(1/n);
end
z = ((1+B)./(1-B))*cc;
x = real(z);
y = imag(z);

function [x,y] = TG( xc, yc, A, B)
N = size(A,2)-1;
zc = complex( xc, yc);
e = complex(zeros(size(zc)),zeros(size(zc)));
for j=1:N+1
    e = e + (A(j)+i*B(j)).*zc.^(1-j);
end
ix = isnan(e);
e(ix) = 0;
z = zc.*exp(e);
x = real(z);
y = imag(z);


function [A,B] = GetTG( N, x, y)
lr = log(sqrt(x.^2 + y.^2));
th = atan2(y,x);
for k=2:size(th,1)
    if (th(k) < th(k-1)) 
        th(k) = th(k) + 2*pi;
    end
end
thi = [th(1:end-1) - 2*pi; th(1:end); th(2:end)+2*pi];
lri = lr([1:end-1,1:end,2:end]);
A = ones(1,N+1);
B = zeros(1,N+1);
Y = complex(zeros(1,2*N), zeros(1,2*N));
tt = 0:pi/N:2*pi;
tt = tt(1:end-1);
Anew = 0*A;
Bnew = B; 
while norm(A-Anew) + norm(B-Bnew) > 1.e-15,
    A = Anew;
    B = Bnew;
    B(1) = th(1) - sum(B(2:N+1));
    B(N+1) = 0;
    Y(1) = 2*N*B(1);
    Y(2:N) = N*(B(2:N)+ i*A(2:N));
    Y(N+1) = 2*N*B(N+1);
    Y(N+2:2*N) = conj(Y(N:-1:2));
    zt = tt + ifft(Y);
    rr = spline(thi, lri, zt);
    Y = fft(rr);
    Anew(1) = real(Y(1))/(2*N);
    Anew(2:N) = real(Y(2:N))/N;
    Anew(N+1) = real(Y(N+1))/(2*N);
    Bnew(1) = B(1);
    Bnew(2:N) = -imag(Y(2:N))/N;
    Bnew(N+1) = 0;
end




