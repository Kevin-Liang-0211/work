% 4th order Adams-Bashforth method 
%   u[n+1] = u[n] + h/24 * (55*f(tn, u[n]) - 59*f(t_n-1, u[n-1]) + 37*f(t_n-2, u[n-2]) - 9*f(t_n-3, u[n-3]))
% 
function [U,t] = adams4(f,u0,t0,tf,N) 

dt=(tf-t0)/N;
t=t0+dt*(0:N);

dim = size(u0);
M = dim(2);

U=zeros(N+1,M);

%U(1,:)=u0;
[U(1:4,:),temp]=rk4(f,u0,t0,t0+3*dt,3); % RK4

fn_3 = f(t(1),U(1,:));
fn_2 = f(t(2),U(2,:));
fn_1 = f(t(3),U(3,:));
for n=4:N
    fn = f(t(n),U(n,:));
    U(n+1,:) = U(n,:) + dt/24*(55*fn-59*fn_1+37*fn_2-9*fn_3);
    fn_3 = fn_2;
    fn_2 = fn_1;
    fn_1 = fn;  
end
