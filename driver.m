clear all;
close all;

np=100;

N=np-1;
h=1/np;
y=(1:1:N)*h;
x=y;
A=matPoisson(np);
f=fPoisson(np);
ex=exact(np);

[X,Y]=meshgrid(x,y);
surf(X,Y,reshapse(f,N,N))