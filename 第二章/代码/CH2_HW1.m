clc; clear; close all;

%% Problem 1: Explicit scheme for diffusion equation
% u_t = a * u_xx,  0 < x < 1, t > 0
% u(x,0) = 1 + sin(pi*x)
% u(0,t) = u(1,t) = 1
% Exact: u(x,t) = 1 + exp(-pi^2*a*t)*sin(pi*x)

a_coeff = 0.1;
xl = 0; xr = 1;
t0 = 0; tf = 1.0;

u0_fun = @(x) 1 + sin(pi*x);
g1 = @(t) 1 + 0*t;
g2 = @(t) 1 + 0*t;
u_exact = @(x, t) 1 + exp(-pi^2*a_coeff*t)*sin(pi*x);

%% Grid parameter sets
dx_list = [0.2,   0.2,  0.2,   0.1,    0.1,     0.1,     0.05,     0.05,      0.05];
dt_list = [0.2,   0.1,  0.05,  0.05,   0.025,   0.0125,  0.0125,   0.00625,   0.003125];
group   = [1,     1,    1,     2,      2,       2,       3,        3,         3];

fprintf('Problem 1: Explicit scheme, a=%.1f, tf=%.1f\n\n', a_coeff, tf);
fprintf('%-6s  %-10s  %-10s  %-8s  %-14s  %-14s\n', ...
    'Group', 'dx', 'dt', 'r', 'L2 error', 'Linf error');
fprintf(repmat('-', 1, 75));
fprintf('\n');

for k = 1:length(dx_list)
    dx = dx_list(k);
    dt = dt_list(k);
    M = round((xr - xl) / dx);
    N = round((tf - t0) / dt);
    r = a_coeff * dt / (dx^2);

    [U, x, t] = diffusion_explicit(u0_fun, g1, g2, a_coeff, xl, xr, t0, tf, M, N);

    ue = zeros(M+1, 1);
    for j = 1:M+1
        ue(j) = u_exact(x(j), tf);
    end

    err = U(:, end) - ue;
    L2_err  = sqrt(sum(err.^2) / (M+1));
    Linf_err = max(abs(err));

    fprintf('(%d)     %-10.4f  %-10.6f  %-8.4f  %-14.6e  %-14.6e\n', ...
        group(k), dx, dt, r, L2_err, Linf_err);
end

%% Plot for group (2) with dx=0.1, dt=0.025 as a representative case
dx = 0.1; dt = 0.025;
M = round((xr - xl) / dx);
N = round((tf - t0) / dt);
[U, x, t] = diffusion_explicit(u0_fun, g1, g2, a_coeff, xl, xr, t0, tf, M, N);

x_fine = linspace(xl, xr, 200);
figure(1)
t_plot = [0, 0.1, 0.5, 1.0];
colors = {'r', 'b', 'g', 'k'};
for p = 1:length(t_plot)
    n_idx = round(t_plot(p) / dt) + 1;
    if n_idx > N+1, n_idx = N+1; end
    plot(x, U(:, n_idx), [colors{p}, '-o'], 'Linewidth', 1.5, 'MarkerSize', 5);
    hold on;
    ue_fine = arrayfun(@(xi) u_exact(xi, t_plot(p)), x_fine);
    plot(x_fine, ue_fine, [colors{p}, '--'], 'Linewidth', 1);
    hold on;
end
legend('t=0 (num)', 't=0 (exact)', 't=0.1 (num)', 't=0.1 (exact)', ...
       't=0.5 (num)', 't=0.5 (exact)', 't=1.0 (num)', 't=1.0 (exact)', ...
       'Location', 'Best')
set(gca, 'FontSize', 12);
title(['Explicit scheme, \Deltax=', num2str(dx), ', \Deltat=', num2str(dt)]);
xlabel('x'); ylabel('u');
print('-dpng', '-r300', 'CH2_HW1.png');
