clc; clear; close all;

%% Problem 3: du/dt = u - t^2 + 1, u(0)=0.5, t in [0,1]
% Exact: u = (t+1)^2 - 0.5*exp(t)
% Methods: 3rd order Kutta, 4th order Runge-Kutta

t0 = 0;
tf = 1.0;
fh = @(t, u) u - t^2 + 1;
u0 = 0.5;
u_exact = @(t) (t+1)^2 - 0.5*exp(t);

%% Part (1): h = 0.2
N = 5;
dt = (tf - t0) / N;

[U_rk3, t] = rk3(fh, u0, t0, tf, N);
[U_rk4, t] = rk4(fh, u0, t0, tf, N);

ue = zeros(N+1, 1);
for n = 1:N+1
    ue(n) = u_exact(t(n));
end

fprintf('Part (1): h = %.1f\n', dt);
fprintf('%-6s  %-14s  %-14s  %-10s  %-14s  %-10s\n', ...
    't', 'Exact', 'RK3', 'Err_RK3', 'RK4', 'Err_RK4');
fprintf(repmat('-', 1, 80));
fprintf('\n');
for n = 1:N+1
    fprintf('%-6.2f  %-14.10f  %-14.10f  %-10.2e  %-14.10f  %-10.2e\n', ...
        t(n), ue(n), U_rk3(n), abs(ue(n)-U_rk3(n)), ...
        U_rk4(n), abs(ue(n)-U_rk4(n)));
end

%% Plot Part (1)
figure(1)
plot(t, U_rk3, '-r*', 'Linewidth', 2);
hold on;
plot(t, U_rk4, '-bo', 'Linewidth', 2);
hold on;
plot(t, ue, 'k', 'Linewidth', 2);
legend('RK3 (Kutta)', 'RK4 (Classical)', 'Exact', 'Location', 'Best')
set(gca, 'FontSize', 14);
title(['h = ', num2str(dt)]);
xlabel('t');
ylabel('u');
print('-dpng', '-r300', 'CH1_HW3_part1.png');

%% Part (2): Convergence order analysis
h_list = [0.2, 0.1, 0.05, 0.025, 0.0125];
err_rk3 = zeros(length(h_list), 1);
err_rk4 = zeros(length(h_list), 1);

for k = 1:length(h_list)
    h = h_list(k);
    N_k = round((tf - t0) / h);
    [U3, t_k] = rk3(fh, u0, t0, tf, N_k);
    [U4, t_k] = rk4(fh, u0, t0, tf, N_k);
    ue_end = u_exact(tf);
    err_rk3(k) = abs(ue_end - U3(end));
    err_rk4(k) = abs(ue_end - U4(end));
end

fprintf('\nPart (2): Convergence order at t = 1\n');
fprintf('%-10s  %-14s  %-10s  %-14s  %-10s\n', ...
    'h', 'Err_RK3', 'Order_RK3', 'Err_RK4', 'Order_RK4');
fprintf(repmat('-', 1, 65));
fprintf('\n');
for k = 1:length(h_list)
    if k == 1
        fprintf('%-10.4f  %-14.6e  %-10s  %-14.6e  %-10s\n', ...
            h_list(k), err_rk3(k), '--', err_rk4(k), '--');
    else
        order_rk3 = log(err_rk3(k-1)/err_rk3(k)) / log(h_list(k-1)/h_list(k));
        order_rk4 = log(err_rk4(k-1)/err_rk4(k)) / log(h_list(k-1)/h_list(k));
        fprintf('%-10.4f  %-14.6e  %-10.4f  %-14.6e  %-10.4f\n', ...
            h_list(k), err_rk3(k), order_rk3, err_rk4(k), order_rk4);
    end
end

%% Plot convergence
figure(2)
loglog(h_list, err_rk3, '-r*', 'Linewidth', 2, 'MarkerSize', 10);
hold on;
loglog(h_list, err_rk4, '-bo', 'Linewidth', 2, 'MarkerSize', 10);
hold on;
loglog(h_list, h_list.^3 * err_rk3(1)/h_list(1)^3, '--r', 'Linewidth', 1);
hold on;
loglog(h_list, h_list.^4 * err_rk4(1)/h_list(1)^4, '--b', 'Linewidth', 1);
legend('RK3 error', 'RK4 error', 'O(h^3) ref', 'O(h^4) ref', 'Location', 'Best')
set(gca, 'FontSize', 14);
xlabel('h');
ylabel('Error at t=1');
title('Convergence Order');
grid on;
print('-dpng', '-r300', 'CH1_HW3_part2.png');
