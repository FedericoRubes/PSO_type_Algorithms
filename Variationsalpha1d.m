clc
clear all
close all

%% Ackley function minimization in dimension 1

B = 0;
C = 0;
f = @(x) -20*exp(-0.2*abs(x-B))-exp(cos(2*pi*(x-B)))+ 20+ exp(1)+ C;
%f = @(x) (x-B).^2 -10*cos(2*pi*(x-B)) +10 +C;


%% Parameters

N      = 100;        % Number of agents
dt     = 1e-1;      % Time step
T      = 10;        % Final time
niter  = T/dt;
eps    = 1e-6;
alphav  = [10,20,30,40,50];        
sigma  = 0.1;
M      = 500;       % Number of samples
lambda = 1;

% Initial position in [-5,5]
X0 = rand(N,1)*10 -5; 

% weight function
waf = @(x) exp(-alpha*f(x));

%% Method
X     = X0; 
wafx  = zeros(N,1);
XM = zeros(1,M);
dec   = zeros(niter,5);

tic
for t=1:M
for k=1:5
    alpha=alphav(k);
    for i=1:niter
    %while norm (X-X(1),2)> eps
        Z =randn(N,1);                     % N N(0,1) random numbers
        wafx = exp(-alpha*f(X));
        m = (X'*wafx)/sum(wafx);
        X = X - dt*lambda*(X-m) + sigma*sqrt(dt)*(abs(X-m).*Z);
        dec(i,k)= dec(i,k)+ (abs(m))^2;
        %dec(i,k)= sqrt(sum(X.^2));
    end
    X    = X0; 
    wafx = zeros(N,1);
end
end
dec=dec/M;
toc

%% plot

% fX0 = f(X0);
% x = linspace(-5, 5, 200);
% fx = f(x);
% 
% figure
% plot (x, fx, '-b', X0, fX0, 'or', XM(1,:), f(XM(1,:)), '*g','linewidth',2)
% xlabel('x')
% ylabel('function value')
% legend('Ackley function','starting points','ending points')
% title('Minimization of Ackley function')


% figure
% plot (x, fx, '-b')
% title('Ackley function')

% figure
% histogram(XM,50)
% axis([-0.01 0.01 0 50])
% title('Histogram')
figure
semilogy (1:niter, dec(:,1),1:niter, dec(:,2),1:niter, dec(:,3),1:niter, dec(:,4),1:niter, dec(:,5),'linewidth',2)
xlabel('iterations')
title('Variations in exponential weight \alpha')
legend('\alpha = 10','\alpha = 20','\alpha = 30','\alpha = 40','\alpha = 50')


