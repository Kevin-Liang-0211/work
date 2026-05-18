clc; clear; close all;

%% Problem 4: 2D diffusion equation with ADI
% u_t = b*(u_xx + u_yy),  0<x<1, 0<y<1, t>0
% u(x,y,0) = sin(pi*x)*sin(pi*y)
% Boundary: all zero
% Exact: u(x,y,t) = exp(-2*b*pi^2*t)*sin(pi*x)*sin(pi*y)
% b = 0.1, compute L2 error at t = 0.4

b = 0.1;
xl = 0; xr = 1;
yb = 0; yt = 1;
t0 = 0; tf = 0.4;

u0_fun = @(x, y) sin(pi*x) * sin(pi*y);
u_exact = @(x, y, t) exp(-2*b*pi^2*t) * sin(pi*x) * sin(pi*y);
f1 = @(y, t) 0; f2 = @(y, t) 0;
g1 = @(x, t) 0; g2 = @(x, t) 0;

dh_list = [0.1, 0.05, 0.01];

fprintf('Problem 4: 2D ADI, b=0.1, L2 error at t=0.4\n\n');
fprintf('%-10s  %-14s\n', 'dx=dy=dt', 'L2 error');
fprintf(repmat('-', 1, 30));
fprintf('\n');

for k = 1:length(dh_list)
    dh = dh_list(k);
    dx = dh; dy = dh; dt = dh;

    nx = round((xr - xl) / dx);
    ny = round((yt - yb) / dy);
    nt = round((tf - t0) / dt);

    x = xl + (0:nx) * dx;
    y = yb + (0:ny) * dy;
    t = t0 + (0:nt) * dt;

    mu_x = b * dt / (dx^2);
    mu_y = b * dt / (dy^2);

    Uold = zeros(nx+1, ny+1);
    Uhalf = zeros(nx+1, ny+1);
    Unew = zeros(nx+1, ny+1);

    for j = 1:ny+1
        for i = 1:nx+1
            Uold(i, j) = u0_fun(x(i), y(j));
        end
    end

    % Thomas algorithm coefficients
    Nx = nx - 1;
    ax = ones(Nx, 1); bx = ones(Nx, 1); cx = ones(Nx, 1);
    bx(:) = 1 + mu_x;
    ax(:) = 0.5 * mu_x;
    cx(:) = 0.5 * mu_x;
    dx_vec = zeros(Nx, 1);

    Ny = ny - 1;
    ay = ones(Ny, 1); by = ones(Ny, 1); cy = ones(Ny, 1);
    by(:) = 1 + mu_y;
    ay(:) = 0.5 * mu_y;
    cy(:) = 0.5 * mu_y;
    dy_vec = zeros(Ny, 1);

    for n = 1:nt
        % Step 1: x-direction sweep
        for j = 1:ny+1
            Uhalf(1, j)    = f1(y(j), t(n+1));
            Uhalf(nx+1, j) = f2(y(j), t(n+1));
        end
        for i = 1:nx+1
            Uhalf(i, 1)    = g1(x(i), t(n+1));
            Uhalf(i, ny+1) = g2(x(i), t(n+1));
        end

        for j = 2:ny
            for i = 2:nx
                dx_vec(i-1) = 0.5*mu_y*Uold(i, j-1) + (1-mu_y)*Uold(i, j) + ...
                              0.5*mu_y*Uold(i, j+1);
            end
            dx_vec(1)    = dx_vec(1)    + 0.5*mu_x*Uhalf(1, j);
            dx_vec(nx-1) = dx_vec(nx-1) + 0.5*mu_x*Uhalf(nx+1, j);
            Uhalf(2:nx, j) = thomas(ax, bx, cx, dx_vec, Nx);
        end

        % Step 2: y-direction sweep
        for j = 1:ny+1
            Unew(1, j)    = f1(y(j), t(n+1));
            Unew(nx+1, j) = f2(y(j), t(n+1));
        end
        for i = 1:nx+1
            Unew(i, 1)    = g1(x(i), t(n+1));
            Unew(i, ny+1) = g2(x(i), t(n+1));
        end

        for i = 2:nx
            for j = 2:ny
                dy_vec(j-1) = 0.5*mu_x*Uhalf(i-1, j) + (1-mu_x)*Uhalf(i, j) + ...
                              0.5*mu_x*Uhalf(i+1, j);
            end
            dy_vec(1)    = dy_vec(1)    + 0.5*mu_y*Uhalf(i, 1);
            dy_vec(ny-1) = dy_vec(ny-1) + 0.5*mu_y*Uhalf(i, ny+1);
            Unew(i, 2:ny) = thomas(ay, by, cy, dy_vec, Ny);
        end

        Uold = Unew;
    end

    % Compute L2 error at tf
    err_sum = 0;
    for j = 1:ny+1
        for i = 1:nx+1
            ue_ij = u_exact(x(i), y(j), tf);
            err_sum = err_sum + (Unew(i, j) - ue_ij)^2;
        end
    end
    L2_err = sqrt(err_sum / ((nx+1)*(ny+1)));

    fprintf('%-10.4f  %-14.6e\n', dh, L2_err);

    % Plot for the finest grid
    if k == length(dh_list)
        Ue = zeros(nx+1, ny+1);
        for j = 1:ny+1
            for i = 1:nx+1
                Ue(i, j) = u_exact(x(i), y(j), tf);
            end
        end

        figure()
        subplot(1, 2, 1);
        surf(x, y, Unew');
        title('ADI numerical solution');
        xlabel('x'); ylabel('y'); zlabel('u');
        set(gca, 'FontSize', 11);

        subplot(1, 2, 2);
        surf(x, y, Ue');
        title('Exact solution');
        xlabel('x'); ylabel('y'); zlabel('u');
        set(gca, 'FontSize', 11);

        sgtitle(['t = ', num2str(tf), ', dx=dy=dt=', num2str(dh)]);
        print('-dpng', '-r300', 'CH2_HW4.png');
    end
end
