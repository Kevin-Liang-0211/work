clc; clear; close all;

%% Problem 5: u' = -5u, u(0) = 1, t in [0,1]
% Exact: u = exp(-5t)
% Methods: 3rd order Adams-Bashforth, 3rd order Adams-Moulton
% Step size h = 0.1, starting steps use RK3

t0 = 0;
tf = 1.0;
N = 10;
dt = (tf - t0) / N;

fh = @(t, u) -5*u;
u0 = 1.0;
u_exact = @(t) exp(-5*t);

%% Compute numerical solutions
[U_ab3, t] = adams3(fh, u0, t0, tf, N);
[U_am3, t] = adams_moulton3(fh, u0, t0, tf, N);

%% Compute exact solution
ue = zeros(N+1, 1);
for n = 1:N+1
    ue(n) = u_exact(t(n));
end

%% Output results
fprintf('Problem 5: u'' = -5u, u(0)=1, h=%.1f\n', dt);
fprintf('%-6s  %-14s  %-14s  %-12s  %-14s  %-12s\n', ...
    't', 'Exact', 'AB3', 'Err_AB3', 'AM3', 'Err_AM3');
fprintf(repmat('-', 1, 80));
fprintf('\n');
for n = 1:N+1
    fprintf('%-6.2f  %-14.10f  %-14.10f  %-12.4e  %-14.10f  %-12.4e\n', ...
        t(n), ue(n), U_ab3(n), abs(ue(n)-U_ab3(n)), ...
        U_am3(n), abs(ue(n)-U_am3(n)));
end

%% Plot
figure(1)
plot(t, U_ab3, '-r*', 'Linewidth', 2);
hold on;
plot(t, U_am3, '-bo', 'Linewidth', 2);
hold on;
plot(t, ue, 'k', 'Linewidth', 2);
legend('Adams-Bashforth 3', 'Adams-Moulton 3', 'Exact', 'Location', 'Best')
set(gca, 'FontSize', 14);
title(['h = ', num2str(dt)]);
xlabel('t');
ylabel('u');
print('-dpng', '-r300', 'CH1_HW5.png');
