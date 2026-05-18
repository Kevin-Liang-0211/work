% Crank-Nicolson scheme for diffusion-reaction equation
% Equation: u_t = D*u_xx + C*u,  xl<=x<=xr, t0<=t<=tf
% u(x,0) = u0(x)
% u(xl,t) = g1(t), u(xr,t) = g2(t)

function [U, x, t] = diffusion_cn_reaction(u0, g1, g2, D, C, xl, xr, t0, tf, nx, nt)

dx = (xr - xl) / nx;
dt = (tf - t0) / nt;
x = xl + (0:nx) * dx;
t = t0 + (0:nt) * dt;

r = D * dt / (2 * dx^2);
c = C * dt / 2;

U = zeros(nx+1, nt+1);
for j = 1:nx+1
    U(j, 1) = u0(x(j));
end
for n = 1:nt+1
    U(1, n) = g1(t(n));
    U(nx+1, n) = g2(t(n));
end

N = nx - 1;
a_vec = ones(N, 1);
b_vec = ones(N, 1);
c_vec = ones(N, 1);
b_vec(:) = 1 + 2*r - c;
a_vec(:) = r;
c_vec(:) = r;
d_vec = zeros(N, 1);

for n = 1:nt
    for j = 2:nx
        d_vec(j-1) = r*U(j-1, n) + (1 - 2*r + c)*U(j, n) + r*U(j+1, n);
    end
    d_vec(1)   = d_vec(1)   + r * U(1, n+1);
    d_vec(N)   = d_vec(N)   + r * U(nx+1, n+1);
    U(2:nx, n+1) = thomas(a_vec, b_vec, c_vec, d_vec, N);
end
