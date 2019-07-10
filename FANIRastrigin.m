clc
clear all
close all

%% Ackley function minimization in dimension d through FA

d = 2;
B = 0;
C = 0;
fr = @(x) sum((x-B).^2 -10*cos(2*pi*(x-B)) +10)/d +C;

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
   I0(h)= fr(X0(h,:));
end

%% Method
X   = X0; 
I   = I0;
Irandom = zeros(N,1);
BestHist = zeros(niter,d);

tic
for k=1:niter 
    perm=randi(N,1,N);
    Xrandom = X(perm,:);
    for h=1:N
        Irandom(h)= fr(Xrandom(h,:));
    end
    r2 = sum((Xrandom-X).^2,2);       
    X = X + (Irandom < I).*( (beta0*exp(-gamma*r2)).*(Xrandom-X)+ exp(-k/niter)*alpha*randn(N,d));
    for h=1:N
        I(h)= fr(X(h,:));
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
    dec(i)= fr(BestHist(i,:));
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



