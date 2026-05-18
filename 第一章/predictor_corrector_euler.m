% Predictor-Corrector Modified Euler (Heun's method)
% Equation: u_t = f(u,t), t0<=t<=tf
% initial condition: u(t0) = u0
%
% Predictor: u_pred = u[n] + h * f(tn, u[n])
% Corrector: u[n+1] = u[n] + h/2 * (f(tn, u[n]) + f(t_{n+1}, u_pred))

function [U, t] = predictor_corrector_euler(f, u0, t0, tf, N)

dt = (tf - t0) / N;
t = t0 + dt * (0:N);

U = zeros(N+1, 1);

U(1) = u0;
for n = 1:N
    fn = f(t(n), U(n));
    u_pred = U(n) + dt * fn;
    U(n+1) = U(n) + 0.5 * dt * (fn + f(t(n+1), u_pred));
end
