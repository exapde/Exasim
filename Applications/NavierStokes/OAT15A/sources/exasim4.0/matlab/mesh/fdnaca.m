function d=fdnaca(p,xd,thick)
%FDNACA Distance function for NACA foil inside a rectangular domain. 
%  D=FDNACA(P,XD,THICK)
%
%      P:         Node positions (Nx2)
%      XD(4):     [x_min, x_max, y_min, y_max] for far field boundary
%      THICK:     Airfoil thickness in percentage
%
%   See also:  NACA, DPOLY, DRECTANGLE, DDIFF
%
th=pi:-pi/200:pi/2;
xt = cos(th)+1; xt(end)=[];  yt=naca(xt,thick);  
xb=fliplr(xt); yb=-naca(xb,thick);
pv=[xt',yt';1,0;xb',yb'];
dfoil=dpoly(p,pv);                          % distance to foil

drec=drectangle(p,xd(1),xd(2),xd(3),xd(4));     % distance to domain edge

d=ddiff(drec,dfoil);

%--------------------------------------------------------------------------

function y=naca(x,thick)
%NACA Equation for a NACA profile. 
%  Y=NACA(X,THICK)
%
%      Y:         y-coordinate
%      X:         x-coordinate
%      THICK:     Airfoil thickness in percentage
%
y=5*0.01*thick*(0.29690*sqrt(x)-0.12600*x-0.35160*x.^2+0.28430*x.^3-0.10150*x.^4);

%--------------------------------------------------------------------------

function h=fhset(p,varargin)
%FHSET Mesh size function for NACA foil. 
%  H=FSET(P)
%
%      P:         Node positions (Nx2)
%      H:         Desired spacing
%
%   See also:  DCIRCLE
%

% --- Positions of the source points ---
s1x=0; s1y=0;       % leading edge
s2x=1; s2y=0;       % trailing edge


% --- Growth parameters ---
h0    = 0.001;       % initial edge lenght
ratio  = 0.02;       % growth ratio
radius = 0.05;         % 0 for a point source

% --- Distance functions ---
h1 = h0 + max(ratio*dcircle(p,s1x,s1y,radius),0);
h2 = h0 + max(ratio*dcircle(p,s2x,s2y,radius),0);
h12 = min(h1,h2);


h4 = h0 + max(ratio*dcircle(p,0.3,0,radius),0);
h5 = h0 + max(ratio*dcircle(p,0.6,0,radius),0);
h45 = min(h4,h5);

h=min(h12,h45);


