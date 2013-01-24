function f=exact(np)
% Builds the right hand side column vector f for (np-1)x(np-1) unknowns
% (interior points) on a grid with grid spacing h=1/np
N=np-1;
h=1/np;
y=(1:1:N)*h;
x=y;
f=zeros(1,N^2);
for j=1:N
    for i=1:N
        f(N*(j-1)+i)=-(x(i)^4-x(i)^2)*(y(j)^4-y(j)^2);
    end
end
f=f';
end