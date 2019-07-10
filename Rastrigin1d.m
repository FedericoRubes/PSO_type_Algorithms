clc
clear all
close all

%% Rastrigin function minimization in dimension 1

B = 0;
C = 0;
f = @(x) (x-B).^2 -10*cos(2*pi*(x-B)) +10 +C;

%% Parameters

N      = 50;        % Number of agents
dt     = 1e-1;      % Time step
T      = 10;        % Final time
niter  = T/dt;
eps    = 1e-6;
alpha  = 50;        
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

for k=1:M
    for i=1:niter
    %while norm (X-X(1),2)> eps
        Z =randn(N,1);                     % N N(0,1) random numbers
        wafx = waf(X);
        m = (X'*wafx)/sum(wafx);
        X = X - dt*lambda*(X-m) + sigma*sqrt(dt)*(norm(X-m,2).*Z);
    end
    XM(:,k)= X(1);
    X0 = rand(N,1)*10 -5;
    X    = X0; 
    wafx = zeros(N,1);
end

%% plot

fX0 = f(X0);
x = linspace(-5, 5, 200);
fx = f(x);

figure
plot (x, fx, '-b', X0, fX0, 'or', XM(1,:), f(XM(1,:)), '*g','linewidth',2)
xlabel('x')
ylabel('function value')
legend('Rastrigin function','starting points','ending points')
title('Minimization of Rastrigin function')

% figure
% plot (x, fx, '-b')
% title('Rastrigin function')

figure
histogram(XM,50)
axis([-0.01 0.01 0 50])
title('Histogram')