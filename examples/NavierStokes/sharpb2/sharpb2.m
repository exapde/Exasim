function [xl, xu] = sharpb2()
% 2D axisymmetric cone profile using polar coordinates for caps
% Nose tip at (0,0), base tip at (1.91,0)

% Parameters (meters)
R  = 0.27;     % base radius
rt = 0.0415;   % nose radius
ax = 0.10;     % base cap length
L  = 1.81;     % cone length (nose tip to base plane)

N = 400;

%% 1) Nose cap (circle) — polar
theta1 = linspace(0, pi/2, N);     % quarter circle
x1 = rt * (1 - cos(theta1));
r1 = rt * sin(theta1);

%% 2) Cone (straight line)
x2 = linspace(rt, L, N);
r2 = rt + (R - rt) * (x2 - rt) / (L - rt);

% %% 3) Base cap (ellipse) — polar-like parameterization
% theta3 = linspace(pi/2, 0, N);
% x3 = L + ax * cos(theta3);
% r3 = R * sin(theta3);

%% Combine upper profile
x = [x1, x2(2:end)];
r = [r1, r2(2:end)];

xl = [x(:) r(:)];

ra = 0.1;   % nose radius
theta1 = linspace(pi, pi/1.5, N);   
x1 = ra * cos(theta1);
r1 = ra * sin(theta1);

x2 = linspace(x1(end), L, N);
r2 = r1(end) + (4*R - r1(end)) * (x2 - x1(end)) / (L - x1(end));

xu = [x1(:) r1(:); x2(:) r2(:)];

figure(1); clf;
hold on;
plot(xl(:,1), xl(:,2), '-r','LineWidth', 2);
plot(xu(:,1), xu(:,2), '-b','LineWidth', 2);
plot(x1(end), r1(end), 'ok');
axis equal;

end