function [x_exp, y_exp] = distribute_points_exponential(x, y, N, lambda_start, lambda_end)
    % DISTRIBUTE_POINTS_EXPONENTIAL interpolates a smooth curve through (x, y)
    % and distributes N points along the curve with exponential clustering at both ends.
    %
    % INPUT:
    %   x - vector of x-coordinates
    %   y - vector of y-coordinates
    %   N - number of points to distribute along the curve
    %   lambda_start - rate parameter for start clustering
    %   lambda_end - rate parameter for end clustering
    %
    % OUTPUT:
    %   x_exp - x-coordinates of exponentially distributed points
    %   y_exp - y-coordinates of exponentially distributed points

    % Ensure x and y are column vectors
    x = x(:);
    y = y(:);
    
    % Interpolate a smooth curve using a spline
    t = cumsum([0; sqrt(diff(x).^2 + diff(y).^2)]);  % Parametric distance (arc length)
    spline_x = spline(t, x);
    spline_y = spline(t, y);

    % Compute the total arc length
    total_length = t(end);

    % Generate exponentially spaced arc length values
    u = linspace(0, 1, N);  % Uniform parameter between 0 and 1

    % Exponential mapping for start and end clustering
    % t_exp_start = total_length * log(1 - u*(1 - exp(-lambda_start))) / -lambda_start;
    t_exp_start = (1 - exp(lambda_start*u))/(1 - exp(lambda_start));
    % t_exp_end = total_length * log(1 - u*(1 - exp(lambda_end))) / lambda_end;
    t_exp_end = (1 - exp(-lambda_end*u))/(1 - exp(-lambda_end));


    % Average the two distributions to get clustering at both ends
    t_exp = cumsum([t_exp_start(1), min(diff(t_exp_start), diff(t_exp_end))]);
    t_exp = total_length * t_exp / t_exp(end);

    % Evaluate the spline at the exponentially spaced arc length values
    x_exp = ppval(spline_x, t_exp);
    y_exp = ppval(spline_y, t_exp);

    % Plot for visualization (optional)
    % figure;
    % plot(x, y, '.-', 'DisplayName', 'Original Points');
    % hold on;
    % plot(x_exp, y_exp, 'ro-', 'DisplayName', ['Exponentially Distributed Points (\lambda_{start} = ', num2str(lambda_start), ', \lambda_{end} = ', num2str(lambda_end), ')']);
    % legend();
    % title(['Exponentially Distributed Points (\lambda_{start} = ', num2str(lambda_start), ', \lambda_{end} = ', num2str(lambda_end), ')']);
    % xlabel('x');
    % ylabel('y');
    % grid on;
    % axis equal;
end


