clc
clear all
close all

%% Ackley function minimization in dimension d through FA

d = 2;
B = 0;
C = 0;
fa = @(x) -20*exp(-(0.2/sqrt(d))*norm(x-B,2))-exp(1/d*sum(cos(2*pi*(x-B))))+ 20+ exp(1)+ C;

%% Parameters

N      = 100;       % Number of agents
niter  = 100;   % between 40 and 100
alpha  = 0.5;        
gamma  = 0.1;   % from 0.01 to 100
beta0  = 1;

% Initial position in [-5,5]
X0 = rand(N,d)*10 -5; 
% Initialize Brightness
I0 = zeros(N,1);
for h=1:N
   I0(h)= fa(X0(h,:));
end

%% Method
X   = X0; 
I   = I0;
IA = zeros(N/2,1);
IB = zeros(N/2,1);
BestHist = zeros(niter,d);

tic
for k=1:niter 
    XA = X(1:N/2,:);
    XB = X(N/2+1:N,:);
    for h=1:N/2
        IA(h)= fa(XA(h,:));
        IB(h)= fa(XB(h,:));
    end
    r2 = sum((XA-XB).^2,2);       
    XA = XA + (IB < IA).*( (beta0*exp(-gamma*r2)).*(XB-XA)+ exp(-k/niter)*alpha*randn(N/2,d));
    XB = XB + (IA < IB).*( (beta0*exp(-gamma*r2)).*(XA-XB)+ exp(-k/niter)*alpha*randn(N/2,d));
    X= [XA;XB];
    X= X(randperm(N),:);
    for h=1:N
        I(h)= fa(X(h,:));
    end
    [m,minpos] = min(I);
    BestHist(k,:) = X(minpos,:);
end
toc

%% Grafici 
t=1:niter;
dec = zeros(niter,1);
dis = zeros(niter,1);
for i=1:niter
    dec(i)= fa(BestHist(i,:));
    dis(i)= norm(BestHist(i,:)-B,1);
end

% Grafico del decadimento
figure
semilogy(t,dec,'linewidth',2.5)
title("Function Value",'FontSize',18)
xlabel('Iterations','FontSize',12)

% Grafico della distanza
figure
semilogy(t,dis,'linewidth',2.5)
title("Error Decay",'FontSize',18)
xlabel('Iterations','FontSize',12)



