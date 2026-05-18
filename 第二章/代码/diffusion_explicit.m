%Explicit scheme for the diffusion equation 
%Equation: u_t = a u_xx, xl<=x<=xr, t0<=t<=tf
%initial condition: u(x,0) = u0(x)
%boundary condition: u(xl,t) = g1(t), u(xr,t) = g2(t) 
function [U,x,t] = diffusion_explicit(u0,g1,g2,a,xl,xr,t0,tf,M,N) 
%input: space interval [xl, xr], time interval [t0, tf],
%       number of space steps M, number of time steps N 
%output: solution [U,x,t]
%Example usage: [U,x,t] = diffusion_explicit(u0,g1,g2,a, 0, 1, 0, 1, 10, 250) 

dx=(xr-xl)/M; dt=(tf-t0)/N; 
x = xl+(0:M)*dx;
t = t0+(0:N)*dt;

mu=dt/(dx*dx);
r=a*mu;
if r > 0.5
    disp('r>0.5, unstable')
end

%initial and boundary condition
U=zeros(M+1,N+1);
for j=1:M+1
    U(j,1)=u0(x(j)); %initial condition
end
for n=1:N+1
    U(1,n)=g1(t(n)); %boundary condition
    U(M+1,n)=g2(t(n));
end

%time stepper
for n=1:N
    for j=2:M
        U(j,n+1)=U(j,n)+r*(U(j-1,n)-2*U(j,n)+U(j+1,n));
    end
end

% U=U';
% 
% % %plot the result
% mesh(x,t,U);
% title('explicit scheme')
% xlabel('space x')
% ylabel('time t')
% zlabel('solution u')
% 
% U=U';

return;
