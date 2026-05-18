% Implicit Euler scheme for ODE 
% Equation: u_t = f(u,t), t0<=t<=tf
% initial condition: u(t0) = u0

% u[n+1] = u[n] + h * f(t_n+1, u[n+1])

function [U,t] = implicit_euler(f,u0,t0,tf,N,eps) 

dt=(tf-t0)/N;
t=t0+dt*(0:N);

U=zeros(N+1,1);

U(1)=u0;
for n=1:N
    U(n+1) = TimeStep(f,t(n),U(n),dt,eps);
end

function unew = TimeStep(f,tn,un,dt,eps)
max_iter = 100; % maximum number for iteration
uold = un; % initial guess
tn1 = tn+dt; % t_{n+1}
for k=1:max_iter
    unew = un+ dt*f(tn1,uold); % update the solution
    if abs(unew-uold)<eps
        break; % converged
    end
    uold = unew;
end
if k>=max_iter
    fprintf('Warning: Divergent!\n');
end

