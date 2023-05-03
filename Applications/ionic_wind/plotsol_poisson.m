function t = plotsol_poisson(shwmsh)
% i is the timestep to plot

% https://www.mathworks.com/matlabcentral/answers/96201-how-do-i-access-a-base-workspace-variable-from-within-a-function
mesh = evalin('base', 'mesh');
sol = evalin('base', 'sol');
sol = sol*-15.15;       % Nondimensionalize and adjust for proper sign

phi = sol(:,1,:);
Er = sol(:,2,:);
Ez = sol(:,3,:);

clf;
t = tiledlayout(2,2);
txt = sprintf('Electrostatic solution');
title(t,txt, 'FontWeight', 'bold')
nexttile
scaplot(mesh,phi,[min(min(phi)) max(max(phi))],0,shwmsh); title('Phi');
nexttile
% scaplot(mesh,Er,[min(min(Er)) max(max(Er))],0,shwmsh); title('Er');
scaplot(mesh,Er,[-1000, 1000],0,shwmsh); title('Er');
nexttile
scaplot(mesh,Ez,[0 2000],0,shwmsh); title('Ez');
nexttile
ndotE = -Ez;    % n=(0,-1)
alpha = 0.5*(tanh(1000000*ndotE)+1);
scaplot(mesh,alpha,[0 1],0,shwmsh); title('alpha');
% xlim([-.005, 0.03]);
% ylim([-.055,-0.03]);
end