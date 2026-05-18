clc; clear; close all;

%% Problem 2: du/dt = u - 2t/u, u(0)=1, t in [0,1]
% Exact: u = sqrt(1 + 2t)
% Methods: Explicit Euler, Modified Euler (iterative), Predictor-Corrector Euler

t0 = 0;
tf = 1.0;
N = 10;
dt = (tf - t0) / N;

fh = @(t, u) u - 2*t/u;
u0 = 1.0;
u_exact = @(t) sqrt(1 + 2*t);
eps = 1.0e-6;

%% Compute numerical solutions
[U1, t] = euler(fh, u0, t0, tf, N);
[U2, t] = modified_euler(fh, u0, t0, tf, N, eps);
[U3, t] = predictor_corrector_euler(fh, u0, t0, tf, N);

%% Compute exact solution
ue = zeros(N+1, 1);
for n = 1:N+1
    ue(n) = u_exact(t(n));
end

%% Plot
figure(1)
plot(t, U1, '-r*', 'Linewidth', 2);
hold on;
plot(t, U2, '-bo', 'Linewidth', 2);
hold on;
plot(t, U3, '-c^', 'Linewidth', 2);
hold on;
plot(t, ue, 'k', 'Linewidth', 2);
legend('Explicit Euler', 'Modified Euler', 'Predictor-Corrector Euler', 'Exact', 'Location', 'Best')
set(gca, 'FontSize', 14);
title(['\Deltat = ', num2str(dt)]);
xlabel('t');
ylabel('u');
print('-dpng', '-r300', 'CH1_HW2.png');

%% Output table
fprintf('Problem 2: du/dt = u - 2t/u, u(0)=1, dt=%.1f\n', dt);
fprintf('%-6s  %-12s  %-12s  %-10s  %-12s  %-10s  %-12s  %-10s\n', ...
    't', 'Exact', 'Euler', 'Err_Euler', 'Mod_Euler', 'Err_Mod', 'PC_Euler', 'Err_PC');
fprintf(repmat('-', 1, 100));
fprintf('\n');
for n = 1:N+1
    fprintf('%-6.2f  %-12.6f  %-12.6f  %-10.2e  %-12.6f  %-10.2e  %-12.6f  %-10.2e\n', ...
        t(n), ue(n), U1(n), abs(ue(n)-U1(n)), ...
        U2(n), abs(ue(n)-U2(n)), ...
        U3(n), abs(ue(n)-U3(n)));
end
