clc
clear all
close all

%% Ackley function minimization in dimension d

d = 2;
B = 0;
C = 0;
f = @(x) -20*exp(-(0.2/sqrt(d))*norm(x-B,2))-exp(1/d*sum(cos(2*pi*(x-B))))+ 20+ exp(1)+ C;

%% Parameters

N      = 100;        % Number of agents
dt     = sqrt(1/N);  % Time step
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
nxi   = zeros(N,1);
for h=1:N
   wafx(h,1)= waf(X(h,:));
end
wafx = wafx/sum(wafx);
normX = zeros(1,niter);

tic
for k=1:niter
    perm=randi(N,1,N);
    Xrandom = X(perm,:);
    X = X + dt*(Xrandom-X).*wafx(perm) + sigma*sqrt(dt)*nxi.*randn(N,d);
    for h=1:N
       wafx(h,1)= waf(X(h,:));
       nxi(h,1) = norm((Xrandom(h,:)-X(h,:)).*wafx(perm(h)),2);
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


