%Crank-Nicolson scheme for the diffusion equation 
%Equation: u_t = a u_xx, xl<=x<=xr, t0<=t<=tf
%initial condition: u(x,0) = u0(x)
%boundary condition: u(xl,t) = g1(t), u(xr,t) = g2(t) 
function [U,x,t] = diffusion_crank_nicolson(u0,g1,g2,a,xl,xr,t0,tf,nx,nt) 
%input: space interval [xl, xr], time interval [t0, tf],
%       number of space steps nx, number of time steps nt
%output: solution [U,x,t]
%Example usage: [U,x,t] = diffusion_crank_nicolson(u0,g1,g2,a,0,1,0,1,10,250) 

dx=(xr-xl)/nx; dt=(tf-t0)/nt; 
x = xl+(0:nx)*dx;
t = t0+(0:nt)*dt;

mu = dt/(dx*dx);
r = a*mu*0.5; 
s = a*mu*0.5;

%initial and boundary condition
U=zeros(nx+1,nt+1);
for j=1:nx+1
    U(j,1)=u0(x(j)); %initial condition
end
for n=1:nt+1
    U(1,n)=g1(t(n)); %boundary condition
    U(nx+1,n)=g2(t(n));
end

% set matrix for the thomas algorithm
N=nx-1;
a=ones(N,1); b=a; c=a; 
b(:)=1+2*r;
a(:)=r;
c(:)=r;
d=zeros(N,1);

%time stepper
for n=1:nt
    for j=2:nx
        %U(j,n+1)=s*U(j-1,n)+(1-2*s)*U(j,n)+s*U(j+1,n); 
        d(j-1)=s*U(j-1,n)+(1-2*s)*U(j,n)+s*U(j+1,n); 
        % note that the spatial index of U is 1 to nx+1,
        % and index of d is 1 to nx-1, which correspond to 2 to nx of U
    end
    d(1)=d(1)+r*U(1,n+1);
    d(nx-1)=d(nx-1)+r*U(nx+1,n+1);
    U(2:nx,n+1)=thomas(a,b,c,d,N);
end

% U=U';
% 
% %plot the result
% mesh(x,t,U);
% title('implicit scheme')
% xlabel('space x')
% ylabel('time t')
% zlabel('solution u')
% 
% U=U';

return;
