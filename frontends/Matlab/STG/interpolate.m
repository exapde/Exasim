function [ltm, Umean, Amean] = interpolate(y, yinput, ltinput, Uinput, Sinput)


% RANS length scale
ltm = interpdata(y, yinput, ltinput);

% % time-average density
% Umean(:,1) = interpdata(y, yinput, Uinput(:,1));
% % time-average x-velocity
% Umean(:,2) = interpdata(y, yinput, Uinput(:,2));
% % time-average y-velocity
% Umean(:,3) = interpdata(y, yinput, Uinput(:,3));
% % time-average z-velocity
% Umean(:,4) = interpdata(y, yinput, Uinput(:,4));
% % time-average temperature
% Umean(:,5) = interpdata(y, yinput, Uinput(:,5));
  
ng = length(y);
nc = size(Uinput,2);
Umean = zeros(ng, nc);
% time-average primitive state
for i = 1:nc
    Umean(:,i) = interpdata(y, yinput, Uinput(:,i));
end

nc = size(Sinput,2);
Amean = zeros(ng, nc);
% time-average Reynolds stress
for i = 1:nc
    Amean(:,i) = interpdata(y, yinput, Sinput(:,i));
end

% Cholesky decomposition of the Reynolds stress tensor
Amean(:,1) = sqrt(Amean(:,1));                               % a11
Amean(:,2) = Amean(:,2)./Amean(:,1);                         % a21
Amean(:,3) = sqrt(Amean(:,3)-Amean(:,2).^2);                 % a22
Amean(:,4) = Amean(:,4)./Amean(:,1);                         % a31
Amean(:,5) = (Amean(:,5)-Amean(:,2).*Amean(:,4))./Amean(:,3); % a32
Amean(:,6) = sqrt(Amean(:,6)-Amean(:,4).^2-Amean(:,5).^2);   % a33


function uoutput = interpdata(y, yinput, uinput)

if isempty(uinput)
    uoutput = 0;
else
    uoutput = interp1(yinput, uinput, y);
end










