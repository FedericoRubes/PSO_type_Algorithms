clc
clear all
close all
%% Time Experiments

%% Ackley function minimization in dimension d

d = 20;
B = 0;
C = 0;
f = @(x) -20*exp(-(0.2/sqrt(d))*norm(x-B,2))-exp(1/d*sum(cos(2*pi*(x-B))))+ 20+ exp(1)+ C;
%f = @(x) sum((x-B).^2 -10*cos(2*pi*(x-B)) +10)/d +C;

%% Parameters

N      = 100;       % Number of agents
dt     = sqrt(1/N);      % Time step
T      = 20;        % Final time
niter  = round(T/dt);
eps    = 1e-3;
alpha  = 50;        
sigma  = 0.1;
lambda = 1;

% Initial position in [-5,5]
X0 = rand(N,d)*10 -5; 

% weight function
waf = @(x) exp(-alpha*f(x));

%% Method
X     = X0; 
Xnew  = X0+1;
wafx  = zeros(N,1);
niter=1;

tic
while norm(Xnew-X,2)>1e-3 && (niter<20000)
    X = Xnew;
    Z =randn(N,d); 
    for h=1:N
       wafx(h,1)= waf(X(h,:));
    end
    wafx = wafx/sum(wafx); 
    m = sum(X.*wafx);
    Xnew = X - dt*lambda*(X-repmat(m,N,1)) + sqrt(dt)* sigma*sqrt(sum((X-repmat(m,N,1)).^2,2)).*Z;
    niter=niter+1;
    
end
toc

%% Grafico se d=2

% for i=1:N
%     fe(i) = f(X0(i,:));
% end
% plot3(X0(:,1),X0(:,2),fe,'ro')
% hold on
% plot3(XM(1,:),XM(2,:),XM(3,:),'*g')
% plot3(0,0,0,'xk')
