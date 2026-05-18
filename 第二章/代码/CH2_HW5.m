clc; clear; close all;

%% Problem 5 (Optional): Advection-Diffusion Equation
% u_t + a*u_x = b*u_xx,  -1 < x < 1,  0 < t < 0.125
% u(x,0) = exp(-10*x^4 / (1 - x^2))
% u(-1,t) = u(1,t) = 0
% FTCS scheme:
%   (U_j^{n+1}-U_j^n)/dt + a*(U_{j+1}^n-U_{j-1}^n)/(2*dx) = b*(U_{j+1}^n-2*U_j^n+U_{j-1}^n)/dx^2
%
% Parameters: a = -10, b = 1

a_adv = -10;
b_diff = 1;
xl = -1; xr = 1;
t0 = 0; tf = 0.125;

u0_fun = @(x) exp(-10*x.^4 ./ (1 - x.^2)) .* (abs(x) < 1);
g1 = @(t) 0;
g2 = @(t) 0;

%% (a) Amplification factor analysis (printed)
fprintf('=== Problem 5: Advection-Diffusion FTCS ===\n\n');
fprintf('(a) Amplification factor:\n');
fprintf('    Assume U_j^n = lambda^n * exp(i*k*j*dx)\n');
fprintf('    lambda = 1 - 2*r*(1 - cos(k*dx)) - i*s*sin(k*dx)\n');
fprintf('    where r = b*dt/dx^2,  s = a*dt/dx\n\n');

%% (b) Stability condition
fprintf('(b) Stability condition: |lambda|^2 <= 1 requires\n');
fprintf('    dt <= dx^2/(2*b)  and  |a|*dx/(2*b) <= 1\n');
fprintf('    Combined: r <= 1/2 and s^2 <= 2r\n\n');

dx = 0.05;
fprintf('    For a=%d, b=%d, dx=%.2f:\n', a_adv, b_diff, dx);
dt_max = dx^2 / (2*b_diff);
fprintf('    dt <= dx^2/(2b) = %.6f\n', dt_max);
cfl_check = abs(a_adv)*dx / (2*b_diff);
fprintf('    |a|*dx/(2b) = %.2f <= 1? %s\n\n', cfl_check, ...
    string(cfl_check <= 1));

%% (c) Numerical computation: dx=0.05, dt=0.001
dx = 0.05;
dt = 0.001;
M = round((xr - xl) / dx);
N = round((tf - t0) / dt);

x = xl + (0:M) * dx;
t_arr = t0 + (0:N) * dt;

r = b_diff * dt / dx^2;
s = a_adv * dt / (2 * dx);

fprintf('Computation: dx=%.4f, dt=%.4f, r=%.4f, s=%.4f\n', dx, dt, r, s);
fprintf('Stability check: r=%.4f <= 0.5? %s\n\n', r, string(r <= 0.5));

U = zeros(M+1, N+1);
for j = 1:M+1
    U(j, 1) = u0_fun(x(j));
end

for n = 1:N
    U(1, n+1) = g1(t_arr(n+1));
    U(M+1, n+1) = g2(t_arr(n+1));
    for j = 2:M
        U(j, n+1) = U(j, n) ...
            - s * (U(j+1, n) - U(j-1, n)) ...
            + r * (U(j+1, n) - 2*U(j, n) + U(j-1, n));
    end
end

%% Evaluate u(x=-0.75, t=0.125)
j_eval = round((-0.75 - xl) / dx) + 1;
u_val = U(j_eval, end);
fprintf('u(x=-0.75, t=0.125) = %.6f\n\n', u_val);

%% Plot at t=0 and t=0.125
figure(1)
plot(x, U(:, 1), 'b-', 'Linewidth', 2);
hold on;
plot(x, U(:, end), 'r-', 'Linewidth', 2);
legend('t = 0', 't = 0.125', 'Location', 'Best');
set(gca, 'FontSize', 14);
xlabel('x'); ylabel('u');
title(['FTCS: a=', num2str(a_adv), ', b=', num2str(b_diff), ...
    ', \Deltax=', num2str(dx), ', \Deltat=', num2str(dt)]);
grid on;
print('-dpng', '-r300', 'CH2_HW5.png');
