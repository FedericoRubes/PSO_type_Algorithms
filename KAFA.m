clc
clear all
close all

%% Ackley function minimization in dimension d through FA

d = 5;
B = 0;
C = 0;
fr = @(x) sum((x-B).^2 -10*cos(2*pi*(x-B)) +10)/d +C;

%% Parameters

N      = 100;   % Number of agents (con 1000 300sec)
niter  = 100;   % between 40 and 100
beta0  = 0.25;
alpha  = sqrt(beta0);        
gamma  = 0.01/d;   % from 0.01 to 100
M=10;
t=1:niter;
dis1 = zeros(niter,1);
dis2 = zeros(niter,1);
dis3 = zeros(niter,1);

% Initial position in [-5,5]
X0 = rand(N,d)*10 -5; 
% Initialize Brightness
I0 = zeros(N,1);
for h=1:N
   I0(h)= fr(X0(h,:));
end

%% Method FA

X1   = X0; 
I1   = I0;
BestHist1 = zeros(niter,d);

tic
for t=1:M
for k=1:niter % Stop condition
    for j=1:N
       r2j = sum((repmat(X1(j,:),N,1)-X1).^2,2);
       %X = X + (I(j)<I).*( (beta0*exp(-gamma*r2j)).*(repmat(X(j,:),N,1)-X)+ alpha*randn(N,d));
       X1 = X1 + (I1(j)<I1).*( (beta0*exp(-gamma*r2j)).*(repmat(X1(j,:),N,1)-X1)+ exp(-k/niter)*alpha*randn(N,d));
       for h=1:N
           I1(h)= fr(X1(h,:));
       end
    end
    [m,minpos] = min(I1);
    BestHist1(k,:) = X1(minpos,:);
    dis1(k)= dis1(k)+norm(BestHist1(k,:)-B,2);
end
end
toc

%% Method FANI
N=5000;
% Initial position in [-5,5]
X0 = rand(N,d)*10 -5; 
% Initialize Brightness
I0 = zeros(N,1);
for h=1:N
   I0(h)= fr(X0(h,:));
end
X2   = X0; 
I2   = I0;
Irandom = zeros(N,1);
BestHist2 = zeros(niter,d);

tic
for t=1:M
for k=1:niter 
    perm=randi(N,1,N);
    Xrandom = X2(perm,:);
    for h=1:N
        Irandom(h)= fr(Xrandom(h,:));
    end
    r2 = sum((Xrandom-X2).^2,2);       
    X2 = X2 + (Irandom < I2).*( (beta0*exp(-gamma*r2)).*(Xrandom-X2)+ exp(-k/niter)*alpha*randn(N,d));
    for h=1:N
        I2(h)= fr(X2(h,:));
    end
    [m,minpos] = min(I2);
    BestHist2(k,:) = X2(minpos,:);
    dis2(k)= dis2(k)+norm(BestHist2(k,:)-B,2);
end
end
toc

%% Method FASN
X3   = X0; 
I3   = I0;
IA = zeros(N/2,1);
IB = zeros(N/2,1);
BestHist3 = zeros(niter,d);

tic
for t=1:M
for k=1:niter 
    XA = X3(1:N/2,:);
    XB = X3(N/2+1:N,:);
    for h=1:N/2
        IA(h)= fr(XA(h,:));
        IB(h)= fr(XB(h,:));
    end
    r2 = sum((XA-XB).^2,2);       
    XA = XA + (IB < IA).*( (beta0*exp(-gamma*r2)).*(XB-XA)+ exp(-k/niter)*alpha*randn(N/2,d));
    XB = XB + (IA < IB).*( (beta0*exp(-gamma*r2)).*(XA-XB)+ exp(-k/niter)*alpha*randn(N/2,d));
    X3= [XA;XB];
    X3= X3(randperm(N),:);
    for h=1:N
        I3(h)= fr(X3(h,:));
    end
    [m,minpos] = min(I3);
    BestHist3(k,:) = X3(minpos,:);
    dis3(k)= dis3(k)+norm(BestHist3(k,:)-B,2);
end
end
toc


%% Grafici 
t=1:niter;
% Grafico della distanza
figure
semilogy(t,dis1/(M*d),'b',t,dis2/(M*d),'r',t,dis3/(M*d),'g','linewidth',2)
title("Mean Error Decay",'FontSize',18)
xlabel('Iterations','FontSize',12)
hold on 

%%
fa = @(x) -20*exp(-(0.2/sqrt(d))*norm(x-B,2))-exp(1/d*sum(cos(2*pi*(x-B))))+ 20+ exp(1)+ C;
dis1 = zeros(niter,1);
dis2 = zeros(niter,1);
dis3 = zeros(niter,1);

N=100;
% Initial position in [-5,5]
X0 = rand(N,d)*10 -5; 
% Initialize Brightness
I0 = zeros(N,1);
for h=1:N
   I0(h)= fr(X0(h,:));
end

%% Method FA
X1   = X0; 
I1   = I0;
BestHist1 = zeros(niter,d);

tic
for t=1:M
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
    dis1(k)= dis1(k)+norm(BestHist1(k,:)-B,2);
end
end
toc

%% Method FANI
N=5000;
% Initial position in [-5,5]
X0 = rand(N,d)*10 -5; 
% Initialize Brightness
I0 = zeros(N,1);
for h=1:N
   I0(h)= fr(X0(h,:));
end
X2   = X0; 
I2   = I0;
Irandom = zeros(N,1);
BestHist2 = zeros(niter,d);

tic
for t=1:M
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
    dis2(k)= dis2(k)+norm(BestHist2(k,:)-B,2);
end
end
toc

%% Method FASN
X3   = X0; 
I3   = I0;
IA = zeros(N/2,1);
IB = zeros(N/2,1);
BestHist3 = zeros(niter,d);

tic
for t=1:M
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
    dis3(k)= dis3(k)+norm(BestHist3(k,:)-B,2);
end
end
toc


%% Grafici 
t=1:niter;
% Grafico della distanza
semilogy(t,dis1/(M*d),'b',t,dis2/(M*d),'r',t,dis3/(M*d),'g','linewidth',2)
title("Error Decay",'FontSize',18)
xlabel('Iterations','FontSize',12)
legend("Rastrigin TI","Rastrigin ANI","Rastrigin ASN","Ackley TI","Ackley ANI","Ackley ASN")







