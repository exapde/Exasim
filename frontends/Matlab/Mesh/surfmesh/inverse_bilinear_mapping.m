function [xi, eta] = inverse_bilinear_mapping(x, y, c, xi0, eta0)
    % Inputs:
    % x, y - Physical coordinates
    % c - Coefficients matrix of size (4x2)
    % Outputs:
    % xi, eta - Reference coordinates

    % Initial guess for (xi, eta)
    xi = xi0;
    eta = eta0;

    % Tolerance and maximum iterations
    tol = 1e-10;
    max_iter = 100;

    for iter = 1:max_iter
        % Evaluate f1 and f2
        f1 = c(1,1) + c(2,1)*xi + c(3,1)*eta + c(4,1)*xi*eta - x;
        f2 = c(1,2) + c(2,2)*xi + c(3,2)*eta + c(4,2)*xi*eta - y;

        % Compute partial derivatives (Jacobian matrix)
        df1_dxi = c(2,1) + c(4,1)*eta;
        df1_deta = c(3,1) + c(4,1)*xi;
        df2_dxi = c(2,2) + c(4,2)*eta;
        df2_deta = c(3,2) + c(4,2)*xi;

        J = [df1_dxi, df1_deta; df2_dxi, df2_deta];

        % Solve for updates (delta_xi, delta_eta)
        F = [f1; f2];
        delta = -J \ F;

        % Update xi and eta
        xi = xi + delta(1);
        eta = eta + delta(2);

        % Check convergence
        if norm(delta) < tol          
            return;
        end
    end

    error('Newton''s method did not converge.');
end
 

