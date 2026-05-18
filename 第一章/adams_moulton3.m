% 3rd order Adams-Moulton method (implicit)
% u[n+1] = u[n] + h/12 * (5*f(t_{n+1}, u[n+1]) + 8*f(t_n, u[n]) - f(t_{n-1}, u[n-1]))
%
% Uses predictor-corrector approach:
%   Predictor: 3rd order Adams-Bashforth
%   Corrector: 3rd order Adams-Moulton
% Starting values computed by RK3

function [U, t] = adams_moulton3(f, u0, t0, tf, N)

dt = (tf - t0) / N;
t = t0 + dt * (0:N);

dim = size(u0);
M = dim(2);

U = zeros(N+1, M);

[U(1:3,:), ~] = rk3(f, u0, t0, t0+2*dt, 2);

fn_2 = f(t(1), U(1,:));
fn_1 = f(t(2), U(2,:));
for n = 3:N
    fn = f(t(n), U(n,:));
    % Predictor (Adams-Bashforth 3)
    U_pred = U(n,:) + dt/12*(23*fn - 16*fn_1 + 5*fn_2);
    % Corrector (Adams-Moulton 3)
    f_pred = f(t(n+1), U_pred);
    U(n+1,:) = U(n,:) + dt/12*(5*f_pred + 8*fn - fn_1);
    fn_2 = fn_1;
    fn_1 = fn;
end
