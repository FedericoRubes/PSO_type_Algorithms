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
    perm=randi(N,1,N);
    Xarandom = Xa(perm,:);
    Xrrandom = Xr(perm,:);
    Xa = Xa + dt*(Xarandom-Xa).*wafax(perm) + sigma*sqrt(dt)*sqrt(sum(((Xarandom-Xa).*wafax(perm)).^2,2)).*randn(N,d);
    Xr = Xr + dt*(Xrrandom-Xr).*wafrx(perm) + sigma*sqrt(dt)*sqrt(sum(((Xrrandom-Xr).*wafrx(perm)).^2,2)).*randn(N,d);
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
plot(t,disa,'b',t,disr,'r')
xlabel('iterations')
ylabel('mean value')
title('Mean error with Asymptotic Nanbu I')
legend('Ackley function','Rastrigin function')
