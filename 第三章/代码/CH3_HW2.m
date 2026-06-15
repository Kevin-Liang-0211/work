clc; clear; close all;

%% Problem 2 (Optional): Burgers Equation - Shock and Rarefaction
% Conservative form:      u_t + (u^2/2)_x = 0
% Non-conservative form:  u_t + u*u_x = 0
%
% Two initial conditions:
%   (1) u0(x) = 2 if x<0, 1 if x>=0   -> shock
%   (2) u0(x) = 1 if x<0, 2 if x>=0   -> centered rarefaction
%
% Exact solutions:
%   (1) u(x,t) = 2 if x < 1.5*t, 1 if x >= 1.5*t          (shock speed 3/2)
%   (2) u(x,t) = 1 if x<t, x/t if t<=x<2t, 2 if x>=2t      (rarefaction fan)
%
% Conservative upwind:     U_j^{n+1} = U_j^n - (dt/dx)*[ (U_j^n)^2/2 - (U_{j-1}^n)^2/2 ]
% Non-conservative upwind: U_j^{n+1} = U_j^n - (dt/dx)*U_j^n*(U_j^n - U_{j-1}^n)
%
% Domain x in [-1,4], dx = 0.01, dt = 0.0025, plot at t = 1.0

xl = -1; xr = 4;
t0 = 0; tf = 1;
dx = 0.01;
dt = 0.0025;

nx = round((xr - xl) / dx);
nt = round((tf - t0) / dt);
x = xl + (0:nx)' * dx;
nu = dt / dx;

fprintf('=== Problem 2: Burgers Equation ===\n');
fprintf('dx = %.3f, dt = %.4f\n\n', dx, dt);

cases = {'Shock', 'Rarefaction'};
for c = 1:2
    if c == 1
        u0 = @(x) 2.0*(x<0) + 1.0*(x>=0);
        uexact = @(x,t) 2.0*(x < 1.5*t) + 1.0*(x >= 1.5*t);
    else
        u0 = @(x) 1.0*(x<0) + 2.0*(x>=0);
        uexact = @(x,t) 1.0*(x<t) + (x/t).*(x>=t & x<2*t) + 2.0*(x>=2*t);
    end

    % initial values
    U_cons = u0(x);
    U_ncon = u0(x);

    % time stepping (wave speed > 0 => left-biased upwind)
    for n = 1:nt
        Uc = U_cons;
        Un = U_ncon;
        for j = 2:nx+1
            % conservative form
            U_cons(j) = Uc(j) - nu*( Uc(j)^2/2 - Uc(j-1)^2/2 );
            % non-conservative form
            U_ncon(j) = Un(j) - nu*Un(j)*( Un(j) - Un(j-1) );
        end
        % left inflow boundary stays fixed (u0(xl))
        U_cons(1) = Uc(1);
        U_ncon(1) = Un(1);
    end

    % exact solution at tf
    Ue = uexact(x, tf);

    figure(c)
    plot(x, Ue,     'k-',  'Linewidth', 2); hold on;
    plot(x, U_cons, 'r--', 'Linewidth', 1.8);
    plot(x, U_ncon, 'b:',  'Linewidth', 2);
    legend('Exact', 'Conservative upwind', 'Non-conservative upwind', ...
        'Location', 'Best');
    axis([xl, xr, 0.5, 2.5]);
    set(gca, 'FontSize', 14);
    xlabel('x'); ylabel('u');
    title([cases{c}, ': t=', num2str(tf), ...
        ', \Deltax=', num2str(dx), ', \Deltat=', num2str(dt)]);
    grid on; box on;
    print('-dpng', '-r300', ['CH3_HW2_case', num2str(c), '.png']);

    % error report
    e_c = U_cons - Ue;  e_n = U_ncon - Ue;
    fprintf('Case %d (%s) at t=1.0:\n', c, cases{c});
    fprintf('  Conservative     L2=%.4e, Linf=%.4e\n', ...
        sqrt(sum(e_c.^2)/(nx+1)), max(abs(e_c)));
    fprintf('  Non-conservative L2=%.4e, Linf=%.4e\n\n', ...
        sqrt(sum(e_n.^2)/(nx+1)), max(abs(e_n)));
end

fprintf('Note: For the shock case, the conservative scheme captures the\n');
fprintf('correct shock speed (3/2), while the non-conservative scheme gives\n');
fprintf('a wrong shock location. For the rarefaction case both work well.\n');
