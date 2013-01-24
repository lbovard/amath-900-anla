function A=matPoisson(np)
% Builds the problem matrix for (np-1)x(np-1) unknowns
% (interior points) on a grid with spacing h=1/np

N=np-1;
%B=[-ones(N^2,1),-ones(N^2,1),4*ones(N^2,1)];
%d=[-N,-1,0];
%F=spdiags(B,d,N^2,N^2);
h=1/np;
A=-sparse(1:(N-1),2:N,1,N,N);
A=A+A';
A=A+sparse(1:N,1:N,4,N,N);
C=A;
for i=1:(N-1)
    C=blkdiag(C,A);
end
D=sparse(1:(N-1)*N,(N+1):N^2,-1,(N)^2,(N)^2);
A=C+D+D';
A=1/h^2*A;


end
