function vnew= GS(A,vold,f)
% perform one Gauss-Seidel iteration on Au=f using initial guess vold

%size of unknown vector
n=length(vold);
%allocate new array
vnew=zeros(1,n)';
for i=1:n
    sum1=0;
    if i~=1
        sum1=sum(A(i,1:(i-1))*vnew(1:(i-1)));
    end
    sum2=sum(A(i,(i+1):n)*vold((i+1):n));
    vnew(i)=(f(i)-sum1-sum2)/A(i,i);
end
end
    