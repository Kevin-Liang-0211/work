%Beam-Warming scheme for the 1d linear advection equation 
%Equation: u_t + a(x,t) u_x = 0, xl<=x<=xr, t0<=t<=tf
%initial condition: u(x,0) = u0(x)
%boundary condition: u(xl,t) = g1(t), or u(xr,t) = g2(t) 
function [U,x,t] = linear_advection_BW(u0,g1,g2,a,xl,xr,t0,tf,nx,nt) 
%input: space interval [xl, xr], time interval [t0, tf],
%       number of space steps nx, number of time steps nt 
%output: solution [U,x,t]
%Example usage: [U,x,t] = linear_advection_upwind(u0,g1,g2,a,0, 1, 0, 1, 10, 250) 

%%%
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
    j=2;
    r = a(x(j),t(n))*nu;
    U(j,n+1)=(1-r)*U(j,n)+r*U(j-1,n);
    for j=3:nx+1
        r = a(x(j),t(n))*nu;
        U(j,n+1)=0.5*(1-r)*(2-r)*U(j,n)+r*(2-r)*U(j-1,n)+...
                 -0.5*r*(1-r)*U(j-2,n);
    end
end

return;
