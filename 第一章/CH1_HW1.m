clc; clear; close all;

%% Problem 1: Euler method for du/dt = u + t
% u(0) = 1, t in [0,1], dt = 0.01
% Exact solution: u = 2*exp(t) - t - 1

t0 = 0;
tf = 1.0;
N = 100;
dt = (tf - t0) / N;

fh = @(t, u) u + t;
u0 = 1.0;
u_exact = @(t) 2*exp(t) - t - 1;

eps = 1.0e-6;

%% Compute numerical solutions
[U1, t] = euler(fh, u0, t0, tf, N);
[U2, t] = implicit_euler(fh, u0, t0, tf, N, eps);
[U3, t] = modified_euler(fh, u0, t0, tf, N, eps);

%% Compute exact solution
ue = zeros(N+1, 1);
for n = 1:N+1
    ue(n) = u_exact(t(n));
end

%% Plot
figure(1)
plot(t, U1, '-r*', 'Linewidth', 1.5, 'MarkerIndices', 1:10:N+1);
hold on;
plot(t, U2, '-bo', 'Linewidth', 1.5, 'MarkerIndices', 1:10:N+1);
hold on;
plot(t, U3, '-c^', 'Linewidth', 1.5, 'MarkerIndices', 1:10:N+1);
hold on;
plot(t, ue, 'k', 'Linewidth', 2);
legend('Explicit Euler', 'Implicit Euler', 'Modified Euler', 'Exact', 'Location', 'Best')
set(gca, 'FontSize', 14);
title(['\Deltat = ', num2str(dt)]);
xlabel('t');
ylabel('u');
print('-dpng', '-r300', 'CH1_HW1.png');

%% Output errors: first 10 and last 10 nodes
fprintf('=== First 10 nodes ===\n');
fprintf('%6s  %12s  %12s  %12s  %12s  %12s  %12s\n', ...
    't', 'Exact', 'Euler', 'Err_Euler', 'Imp_Euler', 'Err_Imp', 'Mod_Euler');
for n = 1:10
    fprintf('%6.4f  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e\n', ...
        t(n), ue(n), U1(n), abs(ue(n)-U1(n)), U2(n), abs(ue(n)-U2(n)), U3(n));
end

fprintf('\n=== Last 10 nodes ===\n');
fprintf('%6s  %12s  %12s  %12s  %12s  %12s  %12s\n', ...
    't', 'Exact', 'Euler', 'Err_Euler', 'Imp_Euler', 'Err_Imp', 'Mod_Euler');
for n = N-8:N+1
    fprintf('%6.4f  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e\n', ...
        t(n), ue(n), U1(n), abs(ue(n)-U1(n)), U2(n), abs(ue(n)-U2(n)), U3(n));
end

%% Errors table (Explicit Euler only as required)
fprintf('\n=== Explicit Euler: First 10 nodes errors ===\n');
fprintf('%8s  %15s  %15s  %12s\n', 't', 'Exact', 'Euler', 'Error');
for n = 1:10
    fprintf('%8.4f  %15.5e  %15.5e  %12.5e\n', t(n), ue(n), U1(n), abs(ue(n)-U1(n)));
end
fprintf('\n=== Explicit Euler: Last 10 nodes errors ===\n');
fprintf('%8s  %15s  %15s  %12s\n', 't', 'Exact', 'Euler', 'Error');
for n = N-8:N+1
    fprintf('%8.4f  %15.5e  %15.5e  %12.5e\n', t(n), ue(n), U1(n), abs(ue(n)-U1(n)));
end
