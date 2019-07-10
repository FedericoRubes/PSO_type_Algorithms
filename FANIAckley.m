clc
%clear all
%close all

%% Ackley function minimization in dimension d through FA

d = 1;
B = 0;
C = 0;
f = @(x) -20*exp(-(0.2/sqrt(d))*norm(x-B,2))-exp(1/d*sum(cos(2*pi*(x-B))))+ 20+ exp(1)+ C;

%% Parameters

N      = 100000;     % Number of agents
niter  = 100;        % between 40 and 100
gamma  = 0.01/(d);   % from 0.01 to 100
beta0  = log10(N)^3/sqrt(N);  % = epsilon
alpha  = sqrt(beta0);        


% Initial position in [-5,5]
X0 = rand(N,d)*10 -5; 
% Initialize Brightness
I0 = zeros(N,1);
for h=1:N
   I0(h)= f(X0(h,:));
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
        Irandom(h)= f(Xrandom(h,:));
    end
    r2 = sum((Xrandom-X).^2,2);       
    X = X + (Irandom < I).*( (beta0*exp(-gamma*r2)).*(Xrandom-X)+ exp(-k/niter)*alpha*randn(N,d));
    for h=1:N
        I(h)= f(X(h,:));
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
    dec(i)= f(BestHist(i,:));
    dis(i)= norm(BestHist(i,:)-B,2)/d;
end

% Grafico del decadimento
% figure
% semilogy(t,dec,'linewidth',2.5)
% title("Function Value",'FontSize',18)
% xlabel('Iterations','FontSize',12)
% 
% % Grafico della distanza
% figure
% semilogy(t,dis,'linewidth',2.5)
% title("Error Decay",'FontSize',18)
% xlabel('Iterations','FontSize',12)




