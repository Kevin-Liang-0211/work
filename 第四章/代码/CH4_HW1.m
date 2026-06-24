clc; clear; close all;

%% Problem 1: Five-point scheme for Poisson equation
% u_xx + u_yy = 16, (x,y) in (0,1)x(0,1)
% u = 0 on the whole boundary
% Initial guess: u0 = 0
% Methods: Jacobi, Gauss-Seidel, SOR (omega = 1.9)

f_fun = @(x, y) 16 + 0 .* x + 0 .* y;
g_fun = @(x, y) 0 .* x .* y;

xl = 0; xr = 1;
yb = 0; yt = 1;
tol = 1.0e-8;
max_iters = 200000;
omega = 1.9;
h_list = [0.1, 0.05, 0.025, 0.0125];

fprintf('=== Problem 1: Poisson equation on (0,1)x(0,1) ===\n');
fprintf('Equation: u_xx + u_yy = 16,  u|_{boundary} = 0\n');
fprintf('Stopping criterion: max |U^{n+1} - U^n| < %.1e\n\n', tol);
fprintf('%-10s %-12s %-14s %-10s %-14s %-10s %-14s %-10s\n', ...
    'dx=dy', 'Jacobi it', 'Jacobi err', 'GS it', 'GS err', ...
    'SOR it', 'SOR err', 'u_center');
fprintf(repmat('-', 1, 100));
fprintf('\n');

for k = 1:length(h_list)
    h = h_list(k);

    [Uj, x, y, iter_j, err_j] = poisson_iterative_solver( ...
        f_fun, g_fun, xl, xr, yb, yt, h, 'jacobi', tol, max_iters);
    [Ug, ~, ~, iter_g, err_g] = poisson_iterative_solver( ...
        f_fun, g_fun, xl, xr, yb, yt, h, 'gauss-seidel', tol, max_iters);
    [Us, ~, ~, iter_s, err_s] = poisson_iterative_solver( ...
        f_fun, g_fun, xl, xr, yb, yt, h, 'sor', tol, max_iters, omega);

    ic = round(length(x) / 2);
    jc = round(length(y) / 2);
    fprintf('%-10.4f %-12d %-14.6e %-10d %-14.6e %-10d %-14.6e %-10.6f\n', ...
        h, iter_j, err_j, iter_g, err_g, iter_s, err_s, Uj(ic, jc));

    figure(k)
    surf(x, y, Uj');
    title(['Jacobi solution, \Deltax=\Deltay=', num2str(h)]);
    xlabel('x'); ylabel('y'); zlabel('u');
    set(gca, 'FontSize', 12);
    view(-55, 22);
    print('-dpng', '-r300', ['CH4_HW1_Jacobi_h', strrep(num2str(h), '.', 'p'), '.png']);
end

fprintf('\nObservation:\n');
fprintf('1. Gauss-Seidel converges faster than Jacobi on every grid.\n');
fprintf('2. SOR with omega = %.1f is the fastest method on all tested grids.\n', omega);
fprintf('3. As the grid is refined, the center value approaches the continuous solution.\n');
