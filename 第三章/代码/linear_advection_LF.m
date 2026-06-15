%Lax-Friedrichs scheme for the 1d linear advection equation
%Equation: u_t + a(x,t) u_x = 0, xl<=x<=xr, t0<=t<=tf
%initial condition: u(x,0) = u0(x)
%boundary condition: u(xl,t) = g1(t), or u(xr,t) = g2(t)
function [U,x,t] = linear_advection_LF(u0,g1,g2,a,xl,xr,t0,tf,nx,nt)
%input: space interval [xl, xr], time interval [t0, tf],
%       number of space steps nx, number of time steps nt
%output: solution [U,x,t]

dx=(xr-xl)/nx; dt=(tf-t0)/nt;
x = xl+(0:nx)*dx;
t = t0+(0:nt)*dt;

nu=dt/dx;

%initial and boundary condition
U=zeros(nx+1,nt+1);
for j=1:nx+1
    U(j,1)=u0(x(j)); %initial condition
end
for n=1:nt+1
    U(1,n)=g1(t(n)); %left boundary condition
    %U(nx+1,n)=g2(t(n)); %right boundary condition
end

%time stepper
for n=1:nt
    for j=2:nx
        r = a(x(j),t(n))*nu;
        U(j,n+1)=0.5*(U(j-1,n)+U(j+1,n))-0.5*r*(U(j+1,n)-U(j-1,n));
    end
    % outflow boundary at right end (first-order upwind)
    r = a(x(nx+1),t(n))*nu;
    U(nx+1,n+1)=(1-r)*U(nx+1,n)+r*U(nx,n);
end

return;
