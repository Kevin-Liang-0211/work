function [U, x, y, iter, max_diff] = poisson_iterative_solver(f_fun, g_fun, xl, xr, yb, yt, h, method, tol, max_iters, omega)
% Solve u_xx + u_yy = f on a rectangle with Dirichlet boundary data
% using Jacobi, Gauss-Seidel, or SOR iteration on the five-point scheme.

if nargin < 11
    omega = 1.9;
end

nx = round((xr - xl) / h);
ny = round((yt - yb) / h);
x = xl + (0:nx) * h;
y = yb + (0:ny) * h;

U = zeros(nx + 1, ny + 1);
for j = 1:ny+1
    U(1, j) = g_fun(x(1), y(j));
    U(nx + 1, j) = g_fun(x(nx + 1), y(j));
end
for i = 1:nx+1
    U(i, 1) = g_fun(x(i), y(1));
    U(i, ny + 1) = g_fun(x(i), y(ny + 1));
end

F = zeros(nx + 1, ny + 1);
for j = 1:ny+1
    for i = 1:nx+1
        F(i, j) = f_fun(x(i), y(j));
    end
end

Unew = U;
max_diff = inf;

for iter = 1:max_iters
    Uold = U;

    switch lower(method)
        case 'jacobi'
            for j = 2:ny
                for i = 2:nx
                    Unew(i, j) = 0.25 * (Uold(i - 1, j) + Uold(i + 1, j) + ...
                        Uold(i, j - 1) + Uold(i, j + 1) - h^2 * F(i, j));
                end
            end
            max_diff = max(max(abs(Unew - Uold)));
            U = Unew;

        case 'gauss-seidel'
            for j = 2:ny
                for i = 2:nx
                    U(i, j) = 0.25 * (U(i - 1, j) + Uold(i + 1, j) + ...
                        U(i, j - 1) + Uold(i, j + 1) - h^2 * F(i, j));
                end
            end
            max_diff = max(max(abs(U - Uold)));

        case 'sor'
            for j = 2:ny
                for i = 2:nx
                    jacobi_val = 0.25 * (U(i - 1, j) + Uold(i + 1, j) + ...
                        U(i, j - 1) + Uold(i, j + 1) - h^2 * F(i, j));
                    U(i, j) = (1 - omega) * Uold(i, j) + omega * jacobi_val;
                end
            end
            max_diff = max(max(abs(U - Uold)));

        otherwise
            error('Unknown method: %s', method);
    end

    if max_diff < tol
        return;
    end
end

warning('%s did not converge within %d iterations.', method, max_iters);
