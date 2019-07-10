clc
clear all
close all

%% Ackley function minimization in dimension d

d = 20;
B = 0;
C = 0;
fa = @(x) -20*exp(-(0.2/sqrt(d))*norm(x-B,2))-exp(1/d*sum(cos(2*pi*(x-B))))+ 20+ exp(1)+ C;
fr = @(x) sum((x-B).^2 -10*cos(2*pi*(x-B)) +10)/d +C;

%% Parameters

N      = 1000;       % Number of agents
dt     = 1e-1;      % Time step
T      = 10;        % Final time
niter  = T/dt;
alpha  = 50;        
sigma  = 0.2;
lambda = 1;
M      = 1;

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
deca   = zeros(niter,1);
decr   = zeros(niter,1);
favet  = zeros(N,1);
frvet  = zeros(N,1);

tic
for t=1:M
    for i=1:niter
        Z =randn(N,d); 
        for h=1:N
           wafax(h,1)= wafa(Xa(h,:));
           wafrx(h,1)= wafr(Xr(h,:));
%            favet(h,1)= norm(ma-B,1);
%            frvet(h,1)= norm(mr-B,1);
        end
        wafax = wafax/sum(wafax); 
        wafrx = wafrx/sum(wafrx); 
        ma = sum(Xa.*wafax);
        mr = sum(Xr.*wafrx);
        Xa = Xa - dt*lambda*(Xa-repmat(ma,N,1)) + sqrt(dt)* sigma*sqrt(sum((Xa-repmat(ma,N,1)).^2,2)).*Z;
        Xr = Xr - dt*lambda*(Xr-repmat(mr,N,1)) + sqrt(dt)* sigma*sqrt(sum((Xr-repmat(mr,N,1)).^2,2)).*Z;
%         deca(i)= sum(favet)/N;
%         decr(i)= sum(frvet)/N;
        deca(i)= deca(i)+norm(ma-B,2)/d;
        decr(i)= decr(i)+norm(mr-B,2)/d;
    end
    X0 = rand(N,d)*10 -5;
    Xa     = X0; 
    Xr     = X0; 
end
toc
deca = deca/M;
decr = decr/M;

%% Grafico andamento media

t=1:1:niter;
figure
semilogy(t,deca,t,decr,'linewidth',2)
legend('Ackley function','Rastrigin function')
xlabel('iterations')
ylabel('mean value')
title('Error decay with total interaction scheme')
