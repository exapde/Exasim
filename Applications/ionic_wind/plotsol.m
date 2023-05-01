function t = plotsol(i, shwmsh)
% i is the timestep to plot

% https://www.mathworks.com/matlabcentral/answers/96201-how-do-i-access-a-base-workspace-variable-from-within-a-function
mesh = evalin('base', 'mesh');
sol = evalin('base', 'sol');

ne = sol(:,1,:,i);
phi = sol(:,2,:,i);
Er = sol(:,4,:,i);
Ez = sol(:,6,:,i);

clf;
t = tiledlayout(3,2);
txt = sprintf('Timestep %d', i);
title(t,txt, 'FontWeight', 'bold')
nexttile
% scaplot(mesh,ne,[min(min(ne)) max(max(ne))],0,shwmsh); title('ne');
scaplot(mesh,ne,[0 0.8],0,shwmsh); title('ne');
nexttile
scaplot(mesh,phi,[min(min(phi)) max(max(phi))],0,shwmsh); title('Phi');
nexttile
% scaplot(mesh,Er,[min(min(Er)) max(max(Er))],0,shwmsh); title('Er');
scaplot(mesh,Er,[-1000, 1000],0,shwmsh); title('Er');
nexttile
% scaplot(mesh,Ez,[min(min(Ez)) max(max(Ez))],0,shwmsh); title('Ez');
scaplot(mesh,Ez,[-1000 1000],0,shwmsh); title('Ez');
nexttile
ndotE = -Ez;    % n=(0,-1)
alpha = 0.5*(tanh(10000*ndotE)+1);
scaplot(mesh,alpha,[-1 1],0,shwmsh); title('alpha');
end