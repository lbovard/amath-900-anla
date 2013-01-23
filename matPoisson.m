function A=matPoisson(np)
% Builds the problem matrix for (np-1)x(np-1) unknowns
% (interior points) on a grid with spacing h=1/np
N=np-1;
A=-sparse(1:(N-1),2:N,1,N,N);
A=A+A';
A=A+sparse(1:N,1:N,4,N,N);
C=A;
for i=1:(N-1)
    C=blkdiag(C,A);
end
D=sparse(1:(N-1)*N,(N+1):N^2,-1,(N)^2,(N)^2);
A=C+D+D';
A=1/np^2*A;
end
