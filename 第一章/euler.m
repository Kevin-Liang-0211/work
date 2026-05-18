% Explicit Euler scheme for ODE 
% Equation: u_t = f(u,t), t0<=t<=tf
% initial condition: u(t0) = u0


% u[n+1] = u[n] + h * f(tn, u[n])

function [U,t] = euler(f,u0,t0,tf,N) 

dt=(tf-t0)/N;
t=t0+dt*(0:N);

U=zeros(N+1,1);

U(1)=u0;
for n=1:N
    U(n+1) = U(n) + dt*f(t(n),U(n));
end
