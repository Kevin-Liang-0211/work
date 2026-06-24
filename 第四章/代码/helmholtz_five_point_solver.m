function [Phi, x, y, max_err] = helmholtz_five_point_solver(lambda, h)
% Solve phi_xx + phi_yy = lambda^2 phi on (0,1)x(0,1)
% with the boundary conditions from the Chapter 4 optional exercise.

mu = sqrt(lambda^2 + pi^2);
phi_exact = @(xv, yv) cos(pi * xv) .* sinh(mu * (1 - yv)) / sinh(mu);
g_left = @(yv) sinh(mu * (1 - yv)) / sinh(mu);
g_right = @(yv) -sinh(mu * (1 - yv)) / sinh(mu);
g_bottom = @(xv) cos(pi * xv);
g_top = @(xv) 0 .* xv;

nx = round(1 / h);
ny = round(1 / h);
x = (0:nx) * h;
y = (0:ny) * h;

Phi = zeros(nx + 1, ny + 1);
for j = 1:ny+1
    Phi(1, j) = g_left(y(j));
    Phi(nx + 1, j) = g_right(y(j));
end
for i = 1:nx+1
    Phi(i, 1) = g_bottom(x(i));
    Phi(i, ny + 1) = g_top(x(i));
end

Nix = nx - 1;
Niy = ny - 1;
N = Nix * Niy;
A = spalloc(N, N, 5 * N);
b = zeros(N, 1);

idx = @(ii, jj) (jj - 1) * Nix + ii;

for j = 1:Niy
    for i = 1:Nix
        row = idx(i, j);
        A(row, row) = -4 - lambda^2 * h^2;

        if i > 1
            A(row, idx(i - 1, j)) = 1;
        else
            b(row) = b(row) - Phi(1, j + 1);
        end

        if i < Nix
            A(row, idx(i + 1, j)) = 1;
        else
            b(row) = b(row) - Phi(nx + 1, j + 1);
        end

        if j > 1
            A(row, idx(i, j - 1)) = 1;
        else
            b(row) = b(row) - Phi(i + 1, 1);
        end

        if j < Niy
            A(row, idx(i, j + 1)) = 1;
        else
            b(row) = b(row) - Phi(i + 1, ny + 1);
        end
    end
end

sol = A \ b;
for j = 1:Niy
    for i = 1:Nix
        Phi(i + 1, j + 1) = sol(idx(i, j));
    end
end

Phi_exact = zeros(nx + 1, ny + 1);
for j = 1:ny+1
    for i = 1:nx+1
        Phi_exact(i, j) = phi_exact(x(i), y(j));
    end
end
max_err = max(max(abs(Phi - Phi_exact)));
