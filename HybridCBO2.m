clc
%clear all
%close all

%% Ackley function minimization in dimension d through FA

d = 5;
B = 0;
C = 0;
%f = @(x) -20*exp(-(0.2/sqrt(d))*norm(x-B,2))-exp(1/d*sum(cos(2*pi*(x-B))))+ 20+ exp(1)+ C;
f = @(x) sum((x-B).^2 -10*cos(2*pi*(x-B)) +10)/d +C;

%% Parameters

N      = 100;       % Number of agents
dt     = 1e-1;      % Time step
T      = 10;        % Final time
niter  = T/dt;

gamma  = 0.01/d;   % from 0.01 to 100
beta0  = 1;
alpha  = 50;        
sigma  = 0.1;
lambda = 1;
alphaf = 0.5;
M      = 3;

t=1:niter;
dish = zeros(niter,1);
dis = zeros(niter,1);
disf = zeros(niter,1);

% Initial position in [-5,5]
X0 = rand(N,d)*10 -5; 

waf = @(x) exp(-alpha*f(x));

% Initialize Brightness
I0 = zeros(N,1);
for h=1:N
   I0(h)= f(X0(h,:));
end

%% Method
for ti=1:M
    
X   = X0; 
I   = I0;
wafx  = zeros(N,1);
m=zeros(niter,d);
index=zeros(N,N);

for k=1:niter % Stop condition
    for h=1:N
        wafx(h,1)= waf(X(h,:));
        index(h,:)=(I(h)>=I).';
    end
    index = index.*wafx';
    wbar = sum(index,2);
    wafx = wafx./wbar;
    for j=1:N
        X = X + (I(j)<I).*( dt*lambda*(wafx(j)).*(repmat(X(j,:),N,1)-X)+...
           sqrt(dt)*exp(-k/niter)* sigma*sqrt(sum((wafx(j).*(repmat(X(j,:),N,1)-X)).^2,2)).*randn(N,d));
    end
    for h=1:N
       %wafx(h,1)= waf(X(h,:));
       I(h)= f(X(h,:));
    end
    m(k,:) = sum(X.*wafx)';
end
mh=m;


% CBO Method
X   = X0; 
wafx  = zeros(N,1);
m=zeros(niter,d);

for i=1:niter
    Z =randn(N,d); 
    for h=1:N
       wafx(h,1)= waf(X(h,:));
    end
    wafx = wafx/sum(wafx); 
    ma = sum(X.*wafx);
    X = X - dt*lambda*(X-repmat(ma,N,1)) + sqrt(dt)* sigma*sqrt(sum((X-repmat(ma,N,1)).^2,2)).*Z;
    m(i,:)=ma;
end
mc=m;

% FA Method
X   = X0; 
I   = I0;
m = zeros(niter,d);

for k=1:niter % Stop condition
    for j=1:N
       r2j = sum((repmat(X(j,:),N,1)-X).^2,2);
       X = X + (I(j)<I).*( (beta0*exp(-gamma*r2j)).*(repmat(X(j,:),N,1)-X)+ exp(-k/niter)*alphaf*randn(N,d));
       for h=1:N
           I(h)= f(X(h,:));
       end
    end
    [mi,minpos] = min(I);
    m(k,:) = X(minpos,:);
end

X0 = rand(N,d)*10 -5; 
for h=1:N
   I0(h)= f(X0(h,:));
end
for i=1:niter
    dish(i)= dish(i)+ norm(mh(i,:)-B,2);
    dis(i)= dis(i)+ norm(mc(i,:)-B,2);
    disf(i)= disf(i)+norm(m(i,:)-B,2);
end
end

%% Grafici 



% Grafico della distanza
figure
semilogy(t,dish/(M*d),t,dis/(M*d),t,disf/(M*d),'linewidth',2)
legend('Hybrid','CBO','FA')
title("Error Decay",'FontSize',18)
xlabel('Iterations','FontSize',12)



