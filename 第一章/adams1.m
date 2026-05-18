% 1st order Adams-Bashforth method */
%  /* Euler method */
%    u[n+1] = u[n] + h * f(tn, u[n])
% 
% 
function [U,t] = adams1(f,u0,t0,tf,N) 

dt=(tf-t0)/N;
t=t0+dt*(0:N);

dim = size(u0);
M = dim(2);

U=zeros(N+1,M);

U(1,:)=u0;
for n=1:N
    U(n+1,:) = U(n,:) + dt*f(t(n),U(n,:));
end
