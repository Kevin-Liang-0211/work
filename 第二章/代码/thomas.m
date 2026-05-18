%Thomas algorithm for solving tridiagonal system
%Equation: Ax=d
% A=[ b1,-c1
%    -a2,b2,-c2
%       -a3,b3,-c3
%        ... ... ...
%             -a_N-1,b_N-1,-c_N-1
%                    -a_N , b_N ]
% d=[d1,d2,d3,...,d_N]
function sol=thomas(a,b,c,d,N) 
%input: a,b,c: entries of matrix A 
%           d: rhs of the equation
%           N: order of the linear system
%output:  sol: the solution
%Example usage: 
% N=5;
% a=ones(N,1); b=a; c=a; d=zeros(N,1);
% b(:)=2; d(1)=1;
% sol=thomas(a,b,c,d,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e=zeros(N,1);
f=zeros(N,1);
x=zeros(N,1);
y=zeros(N,1);

f(1)=b(1);
for i=2:N
    e(i)=-a(i)/f(i-1);
    f(i)=b(i)+e(i)*c(i-1);
end

y(1)=d(1);
for i=2:N
    y(i)=d(i)-e(i)*y(i-1);
end

x(N)=y(N)/f(N);
for i=N-1:-1:1
    x(i)=(y(i)+c(i)*x(i+1))/f(i);
end

sol=x;

return;
