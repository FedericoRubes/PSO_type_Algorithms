clc
clear all
close all

%% Ackley function minimization in dimension d

d = 2;
B = 0;
C = 0;
f = @(x) -20*exp(-(0.2/sqrt(d))*norm(x-B,2))-exp(1/d*sum(cos(2*pi*(x-B))))+ 20+ exp(1)+ C;

%% Parameters

N      = 100;       % Number of agents
dt     = sqrt(1/N);      % Time step
T      = 500;        % Final time
niter  = round(T/dt);
alpha  = 50;        
sigma  = 0.1;
lambda = 1/dt;      % gamma = epsilon = dt

% Initial position in [-5,5]
X0 = rand(N,d)*10 -5; 

% weight function
waf = @(x) exp(-alpha*f(x));

%% Method
X     = X0; 
wafx  = zeros(N,1);
for h=1:N
   wafx(h,1)= waf(X(h,:));
end
wafx = wafx/sum(wafx);
normX = zeros(1,niter);

tic
for k=1:niter
    XA = X(1:N/2,:);
    XB = X(N/2+1:N,:);
    XAnew = XA + dt*(XB-XA).*wafx(N/2+1:N) + sigma*sqrt(dt)*sqrt(sum(((XB-XA).*wafx(N/2+1:N)).^2,2)).*randn(N/2,d);
    XBnew = XB + dt*(XA-XB).*wafx(1:N/2) + sigma*sqrt(dt)*sqrt(sum(((XA-XB).*wafx(1:N/2)).^2,2)).*randn(N/2,d);
    X = [XAnew; XBnew];
    X = X(randperm(N),:);
    for h=1:N
       wafx(h,1)= waf(X(h,:));
    end
    wafx = wafx/sum(wafx);
    normX(k)=norm(X,2);
end
toc

%% Grafico se d=2

figure
plot(X0(:,1),X0(:,2),'ro')
hold on
plot(X(1,:),X(2,:),'*g')
plot(0,0,'xk')

figure
semilogy ([1:niter],normX)
xlabel('iterations')
ylabel('norm of X')
title('Asymptotic Nanbu I')


