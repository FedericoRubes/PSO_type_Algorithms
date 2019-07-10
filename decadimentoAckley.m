clc
clear all
close all

%% Ackley function minimization in dimension d
% Abs to Norm

d = 1;
B = 0;
C = 0;
f = @(x) -20*exp(-(0.2/sqrt(d))*norm(x-B,2))-exp(1/d*sum(cos(2*pi*(x-B))))+ 20+ exp(1)+ C;

%% Parameters

N      = 500;       % Number of agents
dt     = 1e-2;      % Time step
T      = 5;        % Final time
niter  = T/dt;
eps    = 1e-3;
alpha  = 50;        
sigma  = 3;
%M      = 100;       % Number of samples
lambda = 1;

% Initial position in [-5,5]
X0 = rand(N,d)*10 -5; 

% weight function
waf = @(x) exp(-alpha*f(x));

%% Method
X     = X0; 
wafx  = zeros(N,1);
%XM    = zeros(d+1,M);
dec   = zeros(niter,1);
fvet  = zeros(N,1);
index = randi(N,1,5);
exs   = zeros (niter,5);

tic
%for k=1:M
    for i=1:niter
    %while norm (X-repmat(X(1,:),N,1),2)> eps
        Z =randn(N,1);  % NON d
        for h=1:N
           wafx(h,1)= waf(X(h,:));
           fvet(h,1)= f(X(h,:));
        end
        wafx = wafx/sum(wafx); 
        m = sum(X.*wafx);
        X = X - dt*lambda*(X-repmat(m,N,1)) + sqrt(dt)* sigma*(abs(X-repmat(m,N,1)).*Z);
        dec(i)= sum(fvet)/N;
        exs (i,:) = [f(X(index(1),:)),f(X(index(2),:)),f(X(index(3),:)),f(X(index(4),:)),f(X(index(5),:))];
    end
    %XM(1:d,k)= X(1,:)';
    %XM(d+1,k)= f(X(1,:)');
    %X0 = rand(N,d)*10 -5; 
    %X    = X0;
    %wafx = zeros(N,1);
%end
toc

%% Grafico andamento media

t=1:1:niter;
plot(t, dec, 'k--','linewidth',2.5)
hold on
plot(t,exs(:,1),t,exs(:,2),t,exs(:,3),t,exs(:,4),t,exs(:,5))
xlabel('iterations')
ylabel('function value')
title('Rastrigin function decay with iterations')
