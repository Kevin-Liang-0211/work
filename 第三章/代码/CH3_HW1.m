clc; clear; close all;

%% Problem 1: Linear Advection Equation
% u_t + u_x = 0,  0 <= x <= 2,  t > 0
% u(x,0) = u0(x),  u(0,t) = 0
% u0(x) = 1 for 0.4 <= x <= 0.6,  0 otherwise
%
% Solve with Upwind, Lax-Friedrichs, Lax-Wendroff, Beam-Warming schemes.
% dx = 0.01, dt = 0.005  =>  CFL number r = a*dt/dx = 0.5
% Plot solutions at t = 0.5 and t = 1.0
%
% Since a = 1 (constant), the exact solution is u(x,t) = u0(x - t).

%% parameters
xl = 0; xr = 2;
t0 = 0; tf = 1;
dx = 0.01;
dt = 0.005;

nx = round((xr - xl) / dx);
nt = round((tf - t0) / dt);

u0 = @(x) 1.0 * (x >= 0.4 & x <= 0.6);
g1 = @(t) 0*t;
g2 = @(t) 0*t;
a  = @(x,t) 1;

fprintf('=== Problem 1: Linear Advection u_t + u_x = 0 ===\n');
fprintf('dx = %.3f, dt = %.4f, CFL r = a*dt/dx = %.2f\n\n', dx, dt, dt/dx);

%% compute four schemes
[U_up, x, t] = linear_advection_upwind(u0, g1, g2, a, xl, xr, t0, tf, nx, nt);
[U_lf, ~, ~] = linear_advection_LF    (u0, g1, g2, a, xl, xr, t0, tf, nx, nt);
[U_lw, ~, ~] = linear_advection_LW    (u0, g1, g2, a, xl, xr, t0, tf, nx, nt);
[U_bw, ~, ~] = linear_advection_BW    (u0, g1, g2, a, xl, xr, t0, tf, nx, nt);

%% exact solution
Ue = zeros(nx+1, nt+1);
for n = 1:nt+1
    for j = 1:nx+1
        Ue(j,n) = u0(x(j) - t(n));
    end
end

%% plot at selected times
plot_times = [0.5, 1.0];
for k = 1:length(plot_times)
    tk = plot_times(k);
    idx = round((tk - t0)/dt) + 1;

    figure(k)
    plot(x, Ue(:,idx),  'k-',  'Linewidth', 2); hold on;
    plot(x, U_up(:,idx),'r--', 'Linewidth', 1.5);
    plot(x, U_lf(:,idx),'g-.', 'Linewidth', 1.5);
    plot(x, U_lw(:,idx),'b:',  'Linewidth', 2);
    plot(x, U_bw(:,idx),'m-',  'Linewidth', 1.5);
    legend('Exact', 'Upwind', 'Lax-Friedrichs', 'Lax-Wendroff', ...
        'Beam-Warming', 'Location', 'Best');
    axis([0, 2, -0.3, 1.3]);
    set(gca, 'FontSize', 14);
    xlabel('x'); ylabel('u');
    title(['\Deltax=', num2str(dx), ', \Deltat=', num2str(dt), ...
        ', t=', num2str(tk)]);
    grid on; box on;
    print('-dpng', '-r300', ['CH3_HW1_t', num2str(tk*10), '.png']);
end

%% error comparison (L2 and Linf) at t = 1.0
idx = nt + 1;
schemes = {'Upwind', 'Lax-Friedrichs', 'Lax-Wendroff', 'Beam-Warming'};
Us = {U_up, U_lf, U_lw, U_bw};
fprintf('Errors at t = 1.0:\n');
fprintf('%-16s %12s %12s\n', 'Scheme', 'L2 error', 'Linf error');
for s = 1:4
    e = Us{s}(:,idx) - Ue(:,idx);
    L2 = sqrt(sum(e.^2) / (nx+1));
    Linf = max(abs(e));
    fprintf('%-16s %12.4e %12.4e\n', schemes{s}, L2, Linf);
end
