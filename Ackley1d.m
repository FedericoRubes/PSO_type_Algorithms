clc
clear all
close all

%% Ackley function minimization in dimension 1

B = 4;
C = 0;
f = @(x) -20*exp(-0.2*abs(x-B))-exp(cos(2*pi*(x-B)))+ 20+ exp(1)+ C;

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
XMi = zeros(1,M);
XMm = zeros(1,M);
XMf = zeros(1,M);

 wafx = waf(X);
 m = (X'*wafx)/sum(wafx);
        
tic
for k=1:M
    XMi(:,k)= m;
    for i=1:niter
    %while norm (X-X(1),2)> eps
        Z =randn(N,1);                     % N N(0,1) random numbers
        wafx = waf(X);
        m = (X'*wafx)/sum(wafx);
        X = X - dt*lambda*(X-m) + sigma*sqrt(dt)*(abs(X-m).*Z);
        if i==niter/2
            XMm(:,k)= m;
        end
    end
    XMf(:,k)= m;
    X0 = rand(N,1)*10 -5;
    X    = X0; 
    wafx = waf(X);
    m = (X'*wafx)/sum(wafx);
end
toc

%% plot

fX0 = f(X0);
x = linspace(-5, 5, 200);
fx = f(x);

figure
plot (x, fx, '-b', X0, fX0, 'or', XMf(1,:), f(XMf(1,:)), '*g','linewidth',2)
xlabel('x')
ylabel('function value')
legend('Ackley function','starting points','ending points')
title('Minimization of Ackley function')


% figure
% plot (x, fx, '-b')
% title('Ackley function')

%figure
figure
title('m_f^M distribution')
subplot(1,3,1)
histogram(XMi,500)
axis([-1.5 1.5 0 50])
title('t=0')
%figure
subplot(1,3,2)
histogram(XMm,500)
axis([-1.5 1.5 0 50])
title('t=T/2')
%figure
subplot(1,3,3)
histogram(XMf,500)
axis([-1.5 1.5 0 50])
title('t=T')
%axis([-0.01 0.01 0 50])
%title('Histogram')


