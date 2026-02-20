function sharpb2curve()
% Plot the 2D (axisymmetric) x–r profile of the cone geometry from the .geo file.

% Parameters (meters)
cone_length  = 1.81;      % L
base_radius  = 0.27;      % R
tip_radius   = 0.0415;    % rt
base_height  = 0.10;      % ax (semi-axis along x of the base ellipsoid cap)

% Derived
x0 = 0.0;
x1 = cone_length - tip_radius;   % start of spherical nose cap
x2 = cone_length;               % tip location
R  = base_radius;
rt = tip_radius;
ax = base_height;

% Resolution
N1 = 400; N2 = 800; N3 = 200;

%% 1) Base cap: half-ellipse for x in [-ax, 0]
x_base = linspace(-ax, 0, N1);
r_base = R * sqrt(max(0, 1 - (x_base/ax).^2));  % r = R*sqrt(1-(x/ax)^2)

%% 2) Cone frustum: linear radius for x in [0, x1]
x_cone = linspace(0, x1, N2);
r_cone = R + (rt - R) * (x_cone - 0) / (x1 - 0); % linear from R to rt

%% 3) Nose cap: quarter circle for x in [x1, x2]
x_nose = linspace(x1, x2, N3);
r_nose = sqrt(max(0, rt^2 - (x_nose - x1).^2));  % (x-x1)^2 + r^2 = rt^2

%% Assemble outer boundary (upper profile)
x_top = [x_base, x_cone(2:end), x_nose(2:end)];
r_top = [r_base, r_cone(2:end), r_nose(2:end)];

%% Mirror for lower profile (optional, for full cross-section)
x_bot = fliplr(x_top);
r_bot = -fliplr(r_top);

%% Plot
figure; hold on; grid on; box on;
plot(x_top, r_top, 'LineWidth', 2);      % upper boundary
plot(x_bot, r_bot, 'LineWidth', 2);      % lower boundary (mirror)
plot([x_top(1), x_top(end)], [0, 0], 'k--'); % axis line

% Mark key points
plot(-ax, 0, 'ko', 'MarkerFaceColor','k');          % rear tip of base cap
plot(0, R, 'ko', 'MarkerFaceColor','k');            % max radius at x=0
plot(x1, rt, 'ko', 'MarkerFaceColor','k');          % cone-nose junction
plot(x2, 0, 'ko', 'MarkerFaceColor','k');           % nose tip

xlabel('x (m)');
ylabel('r (m)');
title('2D Axisymmetric Geometry: Blunt-base cone with spherical nose');
axis equal;

% Helpful limits
xlim([-(ax*1.1), cone_length*1.02]);
ylim([-(R*1.15), R*1.15]);

legend({'Outer boundary (r\ge0)','Outer boundary (r\le0)','Axis r=0'}, ...
       'Location','best');

end