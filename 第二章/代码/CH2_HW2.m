clc; clear; close all;

%% Problem 2: Implicit and Crank-Nicolson for diffusion equation
% u_t = 4 * u_xx,  0 < x < 1, t > 0
% u(x,0) = exp(-x/2)
% u(0,t) = exp(t),  u(1,t) = exp(t - 1/2)
% Exact: u(x,t) = exp(t - x/2)

a_coeff = 4;
xl = 0; xr = 1;
t0 = 0; tf = 1.0;

u0_fun = @(x) exp(-x/2);
g1 = @(t) exp(t);
g2 = @(t) exp(t - 0.5);
u_exact = @(x, t) exp(t - x/2);

x_eval = 0.5;
t_eval = 1.0;
ue_point = u_exact(x_eval, t_eval);

%% Grid parameters
dh_list = [0.1, 0.05, 0.01];

fprintf('Problem 2: Implicit and C-N, a=4, evaluate at (x,t)=(0.5, 1.0)\n');
fprintf('Exact solution: u(0.5, 1.0) = exp(0.75) = %.10f\n\n', ue_point);

fprintf('%-8s  %-8s  %-16s  %-14s  %-16s  %-14s\n', ...
    'dx', 'dt', 'Implicit', 'Err_Impl', 'Crank-Nicolson', 'Err_CN');
fprintf(repmat('-', 1, 85));
fprintf('\n');

for k = 1:length(dh_list)
    dh = dh_list(k);
    nx = round((xr - xl) / dh);
    nt = round((tf - t0) / dh);

    [U_imp, x, t] = diffusion_implicit(u0_fun, g1, g2, a_coeff, xl, xr, t0, tf, nx, nt);
    [U_cn,  x, t] = diffusion_crank_nicolson(u0_fun, g1, g2, a_coeff, xl, xr, t0, tf, nx, nt);

    j_eval = round((x_eval - xl) / dh) + 1;

    u_imp_val = U_imp(j_eval, end);
    u_cn_val  = U_cn(j_eval, end);

    fprintf('%-8.4f  %-8.4f  %-16.10f  %-14.6e  %-16.10f  %-14.6e\n', ...
        dh, dh, u_imp_val, abs(ue_point - u_imp_val), ...
        u_cn_val, abs(ue_point - u_cn_val));
end

%% Plot for dx=dt=0.05
dh = 0.05;
nx = round((xr - xl) / dh);
nt = round((tf - t0) / dh);
[U_imp, x, t] = diffusion_implicit(u0_fun, g1, g2, a_coeff, xl, xr, t0, tf, nx, nt);
[U_cn,  x, t] = diffusion_crank_nicolson(u0_fun, g1, g2, a_coeff, xl, xr, t0, tf, nx, nt);

x_fine = linspace(xl, xr, 200);

figure(1)
t_vals = [0, 0.25, 0.5, 1.0];
for p = 1:length(t_vals)
    n_idx = round(t_vals(p) / dh) + 1;
    subplot(2, 2, p);
    plot(x, U_imp(:, n_idx), '-b^', 'Linewidth', 1.5);
    hold on;
    plot(x, U_cn(:, n_idx), '-go', 'Linewidth', 1.5);
    hold on;
    ue_fine = arrayfun(@(xi) u_exact(xi, t_vals(p)), x_fine);
    plot(x_fine, ue_fine, 'r-', 'Linewidth', 2);
    legend('Implicit', 'C-N', 'Exact', 'Location', 'Best');
    title(['t = ', num2str(t_vals(p))]);
    xlabel('x'); ylabel('u');
    set(gca, 'FontSize', 11);
end
sgtitle(['dx = dt = ', num2str(dh)]);
print('-dpng', '-r300', 'CH2_HW2.png');
