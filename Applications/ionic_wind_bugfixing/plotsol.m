function t = plotsol(i, shwmsh)
% i is the timestep to plot

% https://www.mathworks.com/matlabcentral/answers/96201-how-do-i-access-a-base-workspace-variable-from-within-a-function
mesh = evalin('base', 'mesh');
sol = evalin('base', 'sol');

clf;
t = tiledlayout(3,2);
txt = sprintf('Timestep %d', i);
title(t,txt, 'FontWeight', 'bold')
nexttile
scaplot(mesh,sol(:,1,:,i),[0 1],0,shwmsh); title('ne');
nexttile
scaplot(mesh,sol(:,2,:,i),[-1 0],0,shwmsh); title('Phi');
nexttile
Er = sol(:,4,:,i);
Ez = sol(:,6,:,i);
scaplot(mesh,Er,[-100 100],0,shwmsh); title('Er');
nexttile
scaplot(mesh,Ez,[-100 100],0,shwmsh); title('Ez');
nexttile
ndotE = 0*Er + -1*Ez;
alpha = 0.5*(tanh(10000*ndotE)+1);
scaplot(mesh,alpha,[-1 1],0,shwmsh); title('alpha');
end