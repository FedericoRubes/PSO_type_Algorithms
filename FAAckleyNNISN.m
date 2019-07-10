clc
clear all
close all

%% Ackley function minimization in dimension d through FA

d = 2;
B = 0;
C = 0;
fa = @(x) -20*exp(-(0.2/sqrt(d))*norm(x-B,2))-exp(1/d*sum(cos(2*pi*(x-B))))+ 20+ exp(1)+ C;

%% Parameters

N      = 100;     % Number of agents
niter  = 1000;   % between 40 and 100
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

%% Method FA
X1   = X0; 
I1   = I0;
BestHist1 = zeros(niter,d);

tic
for k=1:niter % Stop condition
    for j=1:N
       r2j = sum((repmat(X1(j,:),N,1)-X1).^2,2);
       %X = X + (I(j)<I).*( (beta0*exp(-gamma*r2j)).*(repmat(X(j,:),N,1)-X)+ alpha*randn(N,d));
       X1 = X1 + (I1(j)<I1).*( (beta0*exp(-gamma*r2j)).*(repmat(X1(j,:),N,1)-X1)+ exp(-k/niter)*alpha*randn(N,d));
       for h=1:N
           I1(h)= fa(X1(h,:));
       end
    end
    [m,minpos] = min(I1);
    BestHist1(k,:) = X1(minpos,:);
end
toc

%% Method FANI
X2   = X0; 
I2   = I0;
Irandom = zeros(N,1);
BestHist2 = zeros(niter,d);

tic
for k=1:niter 
    perm=randi(N,1,N);
    Xrandom = X2(perm,:);
    for h=1:N
        Irandom(h)= fa(Xrandom(h,:));
    end
    r2 = sum((Xrandom-X2).^2,2);       
    X2 = X2 + (Irandom < I2).*( (beta0*exp(-gamma*r2)).*(Xrandom-X2)+ exp(-k/niter)*alpha*randn(N,d));
    for h=1:N
        I2(h)= fa(X2(h,:));
    end
    [m,minpos] = min(I2);
    BestHist2(k,:) = X2(minpos,:);
end
toc

%% Method FASN
X3   = X0; 
I3   = I0;
IA = zeros(N/2,1);
IB = zeros(N/2,1);
BestHist3 = zeros(niter,d);

tic
for k=1:niter 
    XA = X3(1:N/2,:);
    XB = X3(N/2+1:N,:);
    for h=1:N/2
        IA(h)= fa(XA(h,:));
        IB(h)= fa(XB(h,:));
    end
    r2 = sum((XA-XB).^2,2);       
    XA = XA + (IB < IA).*( (beta0*exp(-gamma*r2)).*(XB-XA)+ exp(-k/niter)*alpha*randn(N/2,d));
    XB = XB + (IA < IB).*( (beta0*exp(-gamma*r2)).*(XA-XB)+ exp(-k/niter)*alpha*randn(N/2,d));
    X3= [XA;XB];
    X3= X3(randperm(N),:);
    for h=1:N
        I3(h)= fa(X3(h,:));
    end
    [m,minpos] = min(I3);
    BestHist3(k,:) = X3(minpos,:);
end
toc


%% Grafici 
t=1:niter;
dec1 = zeros(niter,1);
dec2 = zeros(niter,1);
dec3 = zeros(niter,1);
dis1 = zeros(niter,1);
dis2 = zeros(niter,1);
dis3 = zeros(niter,1);
for i=1:niter
    dec1(i)= fa(BestHist1(i,:));
    dec2(i)= fa(BestHist2(i,:));
    dec3(i)= fa(BestHist3(i,:));
    dis1(i)= norm(BestHist1(i,:)-B,1);
    dis2(i)= norm(BestHist2(i,:)-B,1);
    dis3(i)= norm(BestHist3(i,:)-B,1);
end

% Grafico del decadimento
figure
semilogy(t,dec1,'b',t,dec2,'r',t,dec3,'g','linewidth',2.5)
title("Function Value",'FontSize',18)
xlabel('Iterations','FontSize',12)
legend("Firefly Algorithm","FA with Nanbu I","FA with SN")

% Grafico della distanza
figure
semilogy(t,dis1,'b',t,dis2,'r',t,dis3,'g','linewidth',2.5)
title("Error Decay",'FontSize',18)
xlabel('Iterations','FontSize',12)
legend("Firefly Algorithm","FA with Nanbu I","FA with SN")



