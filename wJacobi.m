function vnew = wJacobi(A,vold,f,omega)
% Perform one weighted Jacobi iteration on Au=f using initial guess vold
% and weight omega

%size of unknown vector
n=length(vold);
%allocate new array
vnew=zeros(1,n)';

%extract main diagonal of A, invert it, store sparsely
Dinv=sparse(diag(1./diag(A)));
%create Jacobi iteration matrix Rj = D^{-1}(L+U)
Rj=Dinv*sparse((-tril(A,-1)-triu(A,1)));
%create weighted Jacobi iteration Rw=(1-w)1+wRj
Rw=(1-omega)*sparse(eye(n,n))+omega*Rj;

%do weighted Jacobi iteration
vnew=Rw*vold+omega*Dinv*f;

end