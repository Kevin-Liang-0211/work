% 3rd order Kutta method
% Equation: u_t = f(u,t), t0<=t<=tf
% initial condition: u(t0) = u0
% 
%    u[n+1] = u[n] + h/6 * (k1 + 4*k2 + k3)
%    k1 = f(tn, u[n])
%    k2 = f(tn+h/2, u[n]+h/2*k1)
%    k3 = f(tn+h, u[n]-h*k1+2h*k2)
% 
function [U,t] = rk3(f,u0,t0,tf,N) 

dt=(tf-t0)/N;
t=t0+dt*(0:N);

dim = size(u0);
M = dim(2);

U=zeros(N+1,M);

U(1,:)=u0;
for n=1:N
    K1 = f(t(n),U(n,:));
    K2 = f(t(n)+0.5*dt,U(n,:)+0.5*dt*K1);
    K3 = f(t(n)+dt,U(n,:)-dt*K1+2*dt*K2);

    U(n+1,:) = U(n,:) + dt/6*(K1+4*K2+K3);
end
