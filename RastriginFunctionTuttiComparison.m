clc
clear all
close all

%% Grafici Ackley function

d = 5;
B = 0;
C = 0;
f = @(x) sum((x-B).^2 -10*cos(2*pi*(x-B)) +10)/d +C;

%% Parameters 

N      = 100;       % Number of agents
T      = 100;        % Final time
alpha  = 50;        
sigma  = 0.1;
M      = 100;       % Number of samples
eps    = 1e-3;

%Total interaction
dtti    = 1e-2;     
niterti = round(T/dtti);
lambda  = 1;

% Binary interaction
dt     = sqrt(1/N);    
niter  = round(T/dt);
gamma = 1/dt;

% Firefly Algorithm
%alpha  = 0.5;        
%gamma  = 0.1;   % from 0.01 to 100
%beta0  = 1;


% Initial position in [-5,5]
X0 = rand(N,d)*10 -5; 
I0 = zeros(N,1);
for h=1:N
   I0(h)= f(X0(h,:));
end

% weight function
waf = @(x) exp(-alpha*f(x));

%% Methods

Xti   = X0; 
Xn1   = X0; 
Xsn   = X0; 
Xfa   = X0;
wafxti  = zeros(N,1);
wafxn1  = zeros(N,1);
wafxsn  = zeros(N,1);
I   = I0;
BestHist = zeros(niter,d);

% Salvo il decadimento della funzione
decti   = zeros(niter,1);
decn1   = zeros(niter,1);
decsn   = zeros(niter,1);
decfa   = zeros(niter,1);

fxti  = zeros(N,1);
fxn1  = zeros(N,1);
fxsn  = zeros(N,1);
fxfa  = zeros(N,1);

% Salvo la distanza dal minimo
disti   = zeros(niter,1);
disn1   = zeros(niter,1);
dissn   = zeros(niter,1);
disfa   = zeros(niter,1);

nti  = zeros(N,1);
nn1  = zeros(N,1);
nsn  = zeros(N,1);
nfa  = zeros(N,1);

tic
    for k=1:niter
        for j=1:N
            r2j = sqrt(sum((repmat(Xfa(j,:),N,1)-Xfa).^2,2));
            Xfa = Xfa + (I(j)<I).*( (exp(-0.1*r2j)).*(repmat(Xfa(j,:),N,1)-Xfa)+ exp(-k/niter)/2*randn(N,d));
            for h=1:N
                I(h)= f(Xfa(h,:));
            end
        end
        [m,minpos] = min(I);
        BestHist(k,:) = Xfa(minpos,:);
        decfa(k)= f(BestHist(k,:));
        disfa(k)= norm(BestHist(k,:)-B,1);
        perm=randi(N,1,N);
        Xn1random = Xn1(perm,:);
        XAa = Xsn(1:N/2,:);
        XBa = Xsn(N/2+1:N,:);
        for h=1:N
           wafxti(h,1)= waf(Xti(h,:));
           wafxn1(h,1)= waf(Xn1(h,:));
           wafxsn(h,1)= waf(Xsn(h,:));
           fxti(h,1)= f(Xti(h,:));
           fxn1(h,1)= f(Xn1(h,:));
           fxsn(h,1)= f(Xsn(h,:));
           nti(h,1)= norm(Xti(h,:)-B,1);
           nn1(h,1)= norm(Xn1(h,:)-B,1);
           nsn(h,1)= norm(Xsn(h,:)-B,1);
        end
        wafxti = wafxti/sum(wafxti);
        wafxn1 = wafxn1/sum(wafxn1);
        wafxsn = wafxsn/sum(wafxsn);
        mti = sum(Xti.*wafxti);
        %Xti = Xti - dtti*lambda*(Xti-repmat(mti,N,1)) + sqrt(dtti)* sigma*sqrt(sum((Xti-repmat(mti,N,1)).^2,2)).*randn(N,d);
        Xti = Xti - dt*lambda*(Xti-repmat(mti,N,1)) + sqrt(dt)* sigma*sqrt(sum((Xti-repmat(mti,N,1)).^2,2)).*randn(N,d);
        Xn1 = Xn1 + dt*(Xn1random-Xn1).*wafxn1(perm) + sigma*sqrt(dt)*sqrt(sum(((Xn1random-Xn1).*wafxn1(perm)).^2,2)).*randn(N,d);
        XAanew = XAa + dt*(XBa-XAa).*wafxsn(N/2+1:N) + sigma*sqrt(dt)*sqrt(sum(((XBa-XAa).*wafxsn(N/2+1:N)).^2,2)).*randn(N/2,d);
        XBanew = XBa + dt*(XAa-XBa).*wafxsn(1:N/2) + sigma*sqrt(dt)*sqrt(sum(((XAa-XBa).*wafxsn(1:N/2)).^2,2)).*randn(N/2,d);
        Xsn = [XAanew; XBanew];
        Xsn = Xsn(randperm(N),:);
        decti(k)  = sum(fxti)/N;
        decn1(k)  = sum(fxn1)/N;
        decsn(k)  = sum(fxsn)/N;
        disti(k)  = sum(nti)/N;
        disn1(k)  = sum(nn1)/N;
        dissn(k)  = sum(nsn)/N;
        end
%     X0 = rand(N,d)*10 -5; 
%     Xti   = X0; 
%     Xn1   = X0; 
%     Xsn   = X0; 
%     wafxti  = zeros(N,1);
%     wafxn1  = zeros(N,1);
%     wafxsn  = zeros(N,1);
toc

%% Grafico
t=1:niter;

% Grafico del decadimento
figure
semilogy(t,decti,'b',t,decn1,'r',t,decsn,'g',t,decfa,'c','linewidth',2.5)
title("Function Mean Value",'FontSize',18)
xlabel('Iterations','FontSize',12)
legend("Total Interaction","Nanbu I","Symmetric Nanbu","Firefly Algorithm")

% Grafico della distanza
figure
semilogy(t,disti,'b',t,disn1,'r',t,dissn,'g',t,disfa,'c','linewidth',2.5)
title("Mean Error Decay",'FontSize',18)
xlabel('Iterations','FontSize',12)
legend("Total Interaction","Nanbu I","Symmetric Nanbu","Firefly Algorithm")


