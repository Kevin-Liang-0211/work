clc; clear; close all;

%% Problem 4: Third-order ODE
% v'''(t) + v''(t) + 4v'(t) + 4v(t) = 4t^2 + 8t - 10
% v(0) = -3, v'(0) = -2, v''(0) = 2
% Exact: v(t) = -sin(2t) + t^2 - 3
% Convert to first-order system: y = [v, v', v'']

t0 = 0;
tf = 1.0;
N = 10;
dt = (tf - t0) / N;

% y1 = v, y2 = v', y3 = v''
% y1' = y2
% y2' = y3
% y3' = 4t^2 + 8t - 10 - 4*y1 - 4*y2 - y3
fh = @(t, y) [y(2), y(3), 4*t^2 + 8*t - 10 - 4*y(1) - 4*y(2) - y(3)];

u0 = [-3, -2, 2];

v_exact = @(t) -sin(2*t) + t^2 - 3;

%% Solve using RK4
[U, t] = rk4(fh, u0, t0, tf, N);

%% Compute exact solution
ue = zeros(N+1, 1);
for n = 1:N+1
    ue(n) = v_exact(t(n));
end

%% Output results
fprintf('Problem 4: v"""+ v"" + 4v'' + 4v = 4t^2+8t-10, dt=%.1f\n', dt);
fprintf('%-6s  %-14s  %-14s  %-12s\n', 't', 'Exact v(t)', 'RK4 v(t)', 'Error');
fprintf(repmat('-', 1, 55));
fprintf('\n');
for n = 1:N+1
    fprintf('%-6.2f  %-14.10f  %-14.10f  %-12.4e\n', ...
        t(n), ue(n), U(n,1), abs(ue(n)-U(n,1)));
end

%% Plot
figure(1)
plot(t, U(:,1), '-bo', 'Linewidth', 2);
hold on;
plot(t, ue, 'k', 'Linewidth', 2);
legend('RK4', 'Exact', 'Location', 'Best')
set(gca, 'FontSize', 14);
title(['v(t), \Deltat = ', num2str(dt)]);
xlabel('t');
ylabel('v');
print('-dpng', '-r300', 'CH1_HW4.png');
