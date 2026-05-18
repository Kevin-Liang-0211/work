clc; clear; close all;

%% Problem 6 (Optional): SIR Model
% dS/dt = -beta * I * S
% dI/dt = beta * I * S - gamma * I
% dR/dt = gamma * I
% N = S + I + R = 10000
% beta = 5.0e-5, gamma = 0.1
% S(0) = 9999, I(0) = 1, R(0) = 0

beta = 5.0e-5;
gamma = 0.1;

fh = @(t, y) [-beta*y(2)*y(1), ...
               beta*y(2)*y(1) - gamma*y(2), ...
               gamma*y(2)];

t0 = 0;
tf = 200;
N = 10000;
dt = (tf - t0) / N;

u0 = [9999, 1, 0];

%% Solve using RK4
[U, t] = rk4(fh, u0, t0, tf, N);

S = U(:, 1);
I_val = U(:, 2);
R = U(:, 3);

%% Find peak of infected
[I_max, idx_max] = max(I_val);
t_peak = t(idx_max);

fprintf('SIR Model Results:\n');
fprintf('Peak infected: I_max = %.2f at t = %.4f\n', I_max, t_peak);
fprintf('At peak: S = %.2f, I = %.2f, R = %.2f\n', S(idx_max), I_val(idx_max), R(idx_max));

%% Plot
figure(1)
plot(t, S, '-b', 'Linewidth', 2);
hold on;
plot(t, I_val, '-r', 'Linewidth', 2);
hold on;
plot(t, R, '-g', 'Linewidth', 2);
hold on;
plot(t_peak, I_max, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
legend('S (Susceptible)', 'I (Infectious)', 'R (Recovered)', ...
    ['Peak I = ', num2str(I_max, '%.1f'), ' at t = ', num2str(t_peak, '%.2f')], ...
    'Location', 'Best')
set(gca, 'FontSize', 14);
xlabel('t');
ylabel('Population');
title('SIR Model');
grid on;
print('-dpng', '-r300', 'CH1_HW6.png');
