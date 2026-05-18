% 2nd order Adams-Bashforth method 
%   u[n+1] = u[n] + h/2 * (3f(tn, u[n])-f(t_n-1, u[n-1]))
% 
function [U,t] = adams2(f,u0,t0,tf,N) 

dt=(tf-t0)/N;
t=t0+dt*(0:N);

dim = size(u0);
M = dim(2);

U=zeros(N+1,M);

U(1,:)=u0;
n=1;
K1 = f(t(n),U(n,:));
K2 = f(t(n)+dt,U(n,:)+dt*K1);
U(2,:)= U(1,:)+0.5*dt*(K1+K2); % RK2

fn_1 = K1;
for n=2:N
    fn = f(t(n),U(n,:));
    U(n+1,:) = U(n,:) + 0.5*dt*(3*fn-fn_1);
    fn_1 = fn;
end
