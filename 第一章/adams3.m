% 3rd order Adams-Bashforth method 
%   u[n+1] = u[n] + h/12 * (23*f(tn, u[n]) - 16*f(t_n-1, u[n-1]) + 5*f(t_n-2, u[n-2]))
% 
function [U,t] = adams3(f,u0,t0,tf,N) 

dt=(tf-t0)/N;
t=t0+dt*(0:N);

dim = size(u0);
M = dim(2);

U=zeros(N+1,M);

%U(1,:)=u0;
[U(1:3,:),temp]=rk3(f,u0,t0,t0+2*dt,2); % RK3

fn_2 = f(t(1),U(1,:));
fn_1 = f(t(2),U(2,:));
for n=3:N
    fn = f(t(n),U(n,:));
    U(n+1,:) = U(n,:) + dt/12*(23*fn-16*fn_1+5*fn_2);
    fn_2 = fn_1;
    fn_1 = fn;
end
