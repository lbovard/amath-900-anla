clear all;
close all;

np=16;
method=1;
omega=4/5;
N=np-1;
h=1/np;
y=(1:1:N)*h;
x=y;
A=matPoisson(np);
f=fPoisson(np);
ex=exact(np);

approx=A\f;
[X,Y]=meshgrid(x,y);
rhs=2*((1-6*X.^2).*(Y.^2).*(1-Y.^2)+(1-6*Y.^2).*(X.^2).*(1-X.^2));

%initial guess
v=rand(1,N^2)';

%compute residual
res_init=norm(f-A*v,2);

if method==0
    disp('Doing GS iteration')
    %do initial GS 
    v=GS(A,v,f);
    res_curr=norm(f-A*v,2);
    res_ratio=res_curr/res_init;
    itrc=1;
    tic
    while(res_ratio>1e-6)
        v=GS(A,v,f);
        res_curr=norm(f-A*v,2);
        res_ratio=res_curr/res_init;
        itrc=itrc+1;
    end
    GStime='GS time took %4.3f seconds\n';
    fprintf(GStime,toc)
    GS_result=norm(v-ex,2);
    GSoutput='GS took %4d iterations with error %1.5f\n';
    fprintf(GSoutput,itrc,GS_result)
    surf(X,Y,reshape(v,N,N),'EdgeColor','none')
    title('GS iteration')
elseif method==1
    disp('Doing Weighted Jacobi')
    %do initial weighted Jacobi
    v=wJacobi(A,v,f,omega);
    res_curr=norm(f-A*v,2);
    res_ratio=res_curr/res_init;
    itrc=1;
    tic
    while(res_ratio>1e-6)
        v=wJacobi(A,v,f,omega);
        res_curr=norm(f-A*v,2);
        res_ratio=res_curr/res_init;
        itrc=itrc+1;
    end
    WJtime='Weighted Jacobi time took %4.3f seconds\n';
    fprintf(WJtime,toc)
    WJ_result=norm(v-ex,2);
    WJoutput='Weighted Jacobi took %4d iterations with error %1.5f\n';
    fprintf(WJoutput,itrc,WJ_result)
    surf(X,Y,reshape(v,N,N),'EdgeColor','none')
    title('Weighted Jacobi iteration')
    
end

figure
%plot the exact solution
surf(X,Y,reshape(ex,N,N),'EdgeColor','none')
title('Exact solution')