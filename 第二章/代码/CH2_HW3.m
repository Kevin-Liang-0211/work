clc; clear; close all;

%% Problem 3: Population growth model
% u_t = D*u_xx + C*u,  0 <= x <= L, t > 0
% u(x,0) = sin^2(pi*x/L)
% u(0,t) = 0,  u(L,t) = 0
% Survival condition: C > pi^2*D/L^2

%% Part (a): L=1, D=1, test C=9.5 and C=10
fprintf('Part (a): L=1, D=1\n');
fprintf('Survival condition: C > pi^2*D/L^2 = %.4f\n\n', pi^2);

L = 1; D_coeff = 1;
xl = 0; xr = L;
dx = 0.01; nx = round(L / dx);
g1 = @(t) 0; g2 = @(t) 0;

C_vals = [9.5, 10];
t_plot = [1, 2, 4, 8];

for ic = 1:length(C_vals)
    C_val = C_vals(ic);
    u0_fun = @(x) sin(pi*x/L).^2;

    figure(ic)
    for p = 1:length(t_plot)
        tf = t_plot(p);
        dt = 0.5 * dx^2 / D_coeff;
        nt = round(tf / dt);
        dt = tf / nt;

        [U, x, t] = diffusion_cn_reaction(u0_fun, g1, g2, D_coeff, C_val, xl, xr, 0, tf, nx, nt);

        subplot(2, 2, p);
        plot(x, U(:, end), 'b-', 'Linewidth', 2);
        title(['t = ', num2str(tf)]);
        xlabel('x'); ylabel('u');
        set(gca, 'FontSize', 11);
        grid on;
    end
    sgtitle(['C = ', num2str(C_val), ' (threshold \approx ', num2str(pi^2, '%.4f'), ')']);
    print('-dpng', '-r300', ['CH2_HW3a_C', num2str(C_val), '.png']);

    fprintf('C = %.1f: ', C_val);
    if C_val > pi^2
        fprintf('C > pi^2, population should SURVIVE\n');
    else
        fprintf('C < pi^2, population should DECAY\n');
    end
end

%% Part (b): C=D=1, survival condition L >= pi
fprintf('\nPart (b): C=D=1\n');
fprintf('Survival condition: L >= pi = %.4f\n\n', pi);

C_val = 1; D_coeff = 1;
L_vals = [3.13, 3.15];
tf = 20;

figure(3)
for il = 1:length(L_vals)
    L = L_vals(il);
    xl = 0; xr = L;
    dx = 0.01; nx = round(L / dx);
    u0_fun = @(x) sin(pi*x/L).^2;
    g1 = @(t) 0; g2 = @(t) 0;

    dt = 0.5 * dx^2 / D_coeff;
    nt = round(tf / dt);
    dt = tf / nt;

    [U, x, t] = diffusion_cn_reaction(u0_fun, g1, g2, D_coeff, C_val, xl, xr, 0, tf, nx, nt);

    max_u = zeros(nt+1, 1);
    for n = 1:nt+1
        max_u(n) = max(U(:, n));
    end

    subplot(1, 2, il);
    plot(t, max_u, 'b-', 'Linewidth', 2);
    xlabel('t'); ylabel('max u(x,t)');
    title(['L = ', num2str(L)]);
    set(gca, 'FontSize', 12);
    grid on;

    fprintf('L = %.2f: ', L);
    if L >= pi
        fprintf('L >= pi, population should SURVIVE\n');
    else
        fprintf('L < pi, population should DECAY\n');
    end
end
sgtitle('Population max density vs time (C=D=1)');
print('-dpng', '-r300', 'CH2_HW3b.png');
