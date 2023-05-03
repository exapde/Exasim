function t = plotsol_alpha(i, shwmsh)
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
% scaplot(mesh,Ez,[min(min(Ez)) max(max(Ez))],0,shwmsh); title('Ez');
scaplot(mesh,Ez,[0 1000],0,shwmsh); title('Ez');
nexttile
ndotE = -Ez;    % n=(0,-1)
alpha1 = 0.5*(tanh(1000000*ndotE)+1);
scaplot(mesh,alpha1,[0 1],0,shwmsh); title('alpha1');
nexttile
alpha2 = 0.5*(tanh(1000000*Ez)+1);
scaplot(mesh,alpha2,[0 1],0,shwmsh); title('alpha2');
end