function y=naca(x,thick)
%NACA Equation for a NACA profile. 
%  Y=NACA(X,THICK)
%
%      Y:         y-coordinate
%      X:         x-coordinate
%      THICK:     Airfoil thickness in percentage
%
y=5*0.01*thick*(0.29690*sqrt(x)-0.12600*x-0.35160*x.^2+0.28430*x.^3-0.10150*x.^4);
