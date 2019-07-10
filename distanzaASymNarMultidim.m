clc
clear all
close all

%% Descrescita della media della funzione f ad ogni iterazione

%% Ackley function minimization in dimension d

d = 2;
B = 0;
C = 0;
fa = @(x) -20*exp(-(0.2/sqrt(d))*norm(x-B,2))-exp(1/d*sum(cos(2*pi*(x-B))))+ 20+ exp(1)+ C;
fr = @(x) sum((x-B).^2 -10*cos(2*pi*(x-B)) +10)/d +C;


%% Parameters

N      = 100;       % Number of agents
dt     = sqrt(1/N);      % Time step
T      = 500;        % Final time
niter  = round(T/dt);
alpha  = 50;        
sigma  = 0.1;
lambda = 1/dt;      % gamma = epsilon = dt

% Initial position in [-5,5]
X0 = rand(N,d)*10 -5; 

% weight function
wafa = @(x) exp(-alpha*fa(x));
wafr = @(x) exp(-alpha*fr(x));

%% Method
Xa     = X0; 
Xr     = X0;
wafax  = zeros(N,1);
wafrx  = zeros(N,1);
for h=1:N
   wafax(h,1)= wafa(Xa(h,:));
   wafrx(h,1)= wafr(Xr(h,:));
end
wafax = wafax/sum(wafax);
wafrx = wafrx/sum(wafrx);

fax  = zeros(N,1);
frx  = zeros(N,1);
Xahist=zeros(niter,2);
Xrhist=zeros(niter,2);
disa   = zeros(niter,1);
disr   = zeros(niter,1);


for k=1:niter
    XAa = Xa(1:N/2,:);
    XBa = Xa(N/2+1:N,:);
    XAanew = XAa + dt*(XBa-XAa).*wafax(N/2+1:N) + sigma*sqrt(dt)*sqrt(sum(((XBa-XAa).*wafax(N/2+1:N)).^2,2)).*randn(N/2,d);
    XBanew = XBa + dt*(XAa-XBa).*wafax(1:N/2) + sigma*sqrt(dt)*sqrt(sum(((XAa-XBa).*wafax(1:N/2)).^2,2)).*randn(N/2,d);
    Xa = [XAanew; XBanew];
    Xa = Xa(randperm(N),:);
    XAr = Xr(1:N/2,:);
    XBr = Xr(N/2+1:N,:);
    XArnew = XAr + dt*(XBr-XAr).*wafrx(N/2+1:N) + sigma*sqrt(dt)*sqrt(sum(((XBr-XAr).*wafrx(N/2+1:N)).^2,2)).*randn(N/2,d);
    XBrnew = XBr + dt*(XAr-XBr).*wafrx(1:N/2) + sigma*sqrt(dt)*sqrt(sum(((XAr-XBr).*wafrx(1:N/2)).^2,2)).*randn(N/2,d);
    Xr = [XArnew; XBrnew];
    Xr = Xr(randperm(N),:);
    for h=1:N
       wafax(h,1)= wafa(Xa(h,:));
       wafrx(h,1)= wafr(Xr(h,:));
       fax(h,1)= norm(Xa(h,:)-B,1);
       frx(h,1)= norm(Xr(h,:)-B,1);
    end
    wafax = wafax/sum(wafax);
    wafrx = wafrx/sum(wafrx);
    Xahist(k,:)=sum(fax)/N;
    Xrhist(k,:)=sum(frx)/N;
    disa(k)= sum(fax)/N;
    disr(k)= sum(frx)/N;
end

%% Graphic
t=1:niter;

figure
plot(t,disa,'b-',t,disr,'r-')
legend('Ackley function','Rastrigin function')
xlabel('iterations')
ylabel('mean value')
title('Mean error with Asymptotic symmetric Nanbu ')

