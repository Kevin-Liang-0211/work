clc; clear; close all;

%% Problem 6 (Optional): Ice Freezing in a Lake
% u_t = alpha * u_xx,  0 < x < 5,  t > 0
% u(x,0) = 4       (initial water temperature 4 C)
% u(0,t) = -15     (lake surface cooled to -15 C)
% du/dx(5,t) = 0   (insulated lake bottom, Neumann BC)
%
% alpha_water = 1.14e-2  m^2/d  (for T >= 0)
% alpha_ice   = 9.98e-2  m^2/d  (for T < 0)
% dx = 0.1 m, dt = 0.1 d

alpha_water = 1.14e-2;
alpha_ice   = 9.98e-2;
xl = 0; xr = 5;
dx = 0.1; dt = 0.1;
nx = round((xr - xl) / dx);

x = xl + (0:nx) * dx;

%% (a) Maximum principle analysis
fprintf('=== Problem 6: Ice Freezing ===\n\n');
fprintf('(a) Maximum principle for C-N requires: alpha*dt/(2*dx^2) <= 1/2\n');
fprintf('    i.e., dt <= dx^2/alpha\n\n');
dt_max_water = dx^2 / alpha_water;
dt_max_ice   = dx^2 / alpha_ice;
fprintf('    For water: dt <= %.4f d\n', dt_max_water);
fprintf('    For ice:   dt <= %.4f d\n', dt_max_ice);
fprintf('    dt = %.1f d: water OK=%s, ice OK=%s\n\n', dt, ...
    string(dt <= dt_max_water), string(dt <= dt_max_ice));

s_water = alpha_water * dt / (2 * dx^2);
s_ice   = alpha_ice   * dt / (2 * dx^2);
fprintf('    s_water = %.4f, s_ice = %.4f\n', s_water, s_ice);
fprintf('    1-2*s_water = %.4f >= 0? %s\n', 1-2*s_water, string(1-2*s_water >= 0));
fprintf('    1-2*s_ice   = %.4f >= 0? %s\n\n', 1-2*s_ice, string(1-2*s_ice >= 0));

%% (b) Compute with variable alpha until ice reaches 1m
% Ice region: T < 0 => alpha = alpha_ice
% Water region: T >= 0 => alpha = alpha_water

U = 4 * ones(nx+1, 1);
U(1) = -15;

max_days = 5000;
ice_depth = 0;
t_freeze = 0;

for n = 1:round(max_days / dt)
    tn = n * dt;

    alpha_vec = zeros(nx+1, 1);
    for j = 1:nx+1
        if U(j) < 0
            alpha_vec(j) = alpha_ice;
        else
            alpha_vec(j) = alpha_water;
        end
    end

    Unew = U;

    % Interior points: C-N with local alpha
    % Use average alpha at each point for the time step
    N_sys = nx - 1;
    a_lo = zeros(N_sys, 1);
    b_diag = zeros(N_sys, 1);
    c_up = zeros(N_sys, 1);
    d_rhs = zeros(N_sys, 1);

    for j = 2:nx
        idx = j - 1;
        alpha_j = alpha_vec(j);
        s = alpha_j * dt / (2 * dx^2);

        b_diag(idx) = 1 + 2*s;
        a_lo(idx) = s;
        c_up(idx) = s;
        d_rhs(idx) = s*U(j-1) + (1 - 2*s)*U(j) + s*U(j+1);
    end

    % Left BC: Dirichlet u(0,t) = -15
    d_rhs(1) = d_rhs(1) + a_lo(1) * (-15);

    % Right BC: Neumann du/dx(5,t) = 0
    % Ghost point: U(nx+2) = U(nx)  =>  at j=nx: s*U(nx+1) term uses U(nx-1)
    % Modify last equation
    j = nx;
    idx = j - 1;
    alpha_j = alpha_vec(j);
    s = alpha_j * dt / (2 * dx^2);
    b_diag(idx) = 1 + 2*s;
    a_lo(idx) = s;
    c_up(idx) = 0;
    d_rhs(idx) = s*U(j-1) + (1 - 2*s)*U(j) + s*U(j+1);
    % For Neumann: U(nx+1)^{n+1} = U(nx-1)^{n+1}, so add s to a_lo and adjust
    % Actually, modify: at j=nx, the stencil uses j+1 = nx+1 (ghost) = j-1 = nx-1
    % LHS: -s*U_{nx-1}^{n+1} + (1+2s)*U_{nx}^{n+1} - s*U_{nx+1}^{n+1}
    %     = -s*U_{nx-1}^{n+1} + (1+2s)*U_{nx}^{n+1} - s*U_{nx-1}^{n+1}
    %     = -2s*U_{nx-1}^{n+1} + (1+2s)*U_{nx}^{n+1}
    % RHS: s*U_{nx-1}^n + (1-2s)*U_{nx}^n + s*U_{nx+1}^n
    %     = s*U_{nx-1}^n + (1-2s)*U_{nx}^n + s*U_{nx-1}^n
    %     = 2s*U_{nx-1}^n + (1-2s)*U_{nx}^n
    a_lo(idx) = 2*s;
    c_up(idx) = 0;
    d_rhs(idx) = 2*s*U(j-1) + (1 - 2*s)*U(j);

    % Also handle the point j = nx+1 (boundary): apply Neumann after solve
    sol = thomas(a_lo, b_diag, c_up, d_rhs, N_sys);
    Unew(2:nx) = sol;
    Unew(1) = -15;
    Unew(nx+1) = Unew(nx);

    U = Unew;

    % Check ice depth: find deepest x where T < 0
    ice_idx = find(U < 0);
    if ~isempty(ice_idx)
        current_depth = x(ice_idx(end));
    else
        current_depth = 0;
    end

    if current_depth >= 1.0 && t_freeze == 0
        t_freeze = tn;
        U_freeze = U;
        fprintf('Ice reaches 1m at t = %.1f days\n', t_freeze);
        break;
    end
end

if t_freeze == 0
    fprintf('Ice did not reach 1m within %d days\n', max_days);
    t_freeze = n * dt;
    U_freeze = U;
end

%% Plot temperature distribution at freezing time
figure(1)
plot(x, U_freeze, 'b-', 'Linewidth', 2);
hold on;
plot([xl, xr], [0, 0], 'k--', 'Linewidth', 1);
hold on;
% Mark the ice-water interface (find zero-crossing)
ice_boundary = NaN;
for jj = 1:nx
    if U_freeze(jj) < 0 && U_freeze(jj+1) >= 0
        frac = (0 - U_freeze(jj)) / (U_freeze(jj+1) - U_freeze(jj));
        ice_boundary = x(jj) + frac * dx;
        break;
    end
end
if ~isnan(ice_boundary)
    plot(ice_boundary, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
end
set(gca, 'FontSize', 14);
xlabel('x (m)');
ylabel('Temperature (°C)');
title(['Temperature at t = ', num2str(t_freeze, '%.1f'), ' days']);
legend('Temperature', '0°C line', 'Ice-water interface', 'Location', 'Best');
grid on;
print('-dpng', '-r300', 'CH2_HW6.png');

fprintf('\nFreezing 1m of ice takes approximately %.1f days\n', t_freeze);
