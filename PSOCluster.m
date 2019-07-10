clc
clear all
close all

%% Parameters

Nd = 2;    % Dimension of the space
No = 300;   % Number of objects
Nc = 5;    % number of clusters
Np = 50;   % Number of particles in PCO
tmax = 100;

X0 = 2*randn (Nc, Nd, Np);  % contains centroids
O = zeros(No, Nd);  % contains the objects
for i=0:No/Nc-1
   O(Nc*i+3,:)= randn(1,Nd)+[4,4] ;
   O(Nc*i+1,:)= randn(1,Nd)+[-4,-4];
   O(Nc*i+2,:)= randn(1,Nd);
   O(Nc*i+5,:)= randn(1,Nd)+[4,-4];
   O(Nc*i+4,:)= randn(1,Nd)+[-4,4];
end

c0 = 0.72;
c1 = 1.49;
c2 = 1.49;
V= zeros (Nc, Nd, Np);


%% Method
dist = zeros(1, Nc);  % distance of 
card = zeros(Nc,Np);
valf = zeros(Nc,Np);

J = zeros (Np+1, 1);        % fitness for every particle + best
Jold = zeros (Np+1, 1)+10^5;        % fitness for every particle + best
actualbest = zeros (Nc,Nd);       % contain global (1) best
lbest = zeros (Nc, Nd, Np); % contain local(Np) best
bestval = zeros (1,Np+1);   % contain local(Np) + global (1) best values
njk = zeros (Nc, Nd, Np);
opt = zeros(tmax,1);

X=X0;
% Insert Kmeans
[idx,C,sumd]=kmeans(O,Nc);
X(:,:,1)= C;

for t=1:tmax
    for i=1:Np
        for j=1:No
            for k=1:Nc
                % Distanza Zj, Mij
                dist(k) = norm(O(j,:)-X(k,:,i),2);                
            end
            [mmin, mminpos] = min(dist);
            card (mminpos,i)= card (mminpos,i)+1;
            valf (mminpos,i)= valf (mminpos,i)+mmin;
        end
        J(i) = sum(valf(:,i)./(card(:,i)+(card(:,i)==0)))/Nc;
    end
    card = zeros(Nc,Np);
    valf = zeros(Nc,Np);
    [mmmin,mmminpos] = min(J(1:Np));
    opt(t)=mmminpos;
    J(Np+1) = mmmin;
    % Change L/G Best
    if J(Np+1)<Jold(Np+1)
        actualbest = X(:,:,mmminpos);
    end
    asd= J(1:Np) < Jold(1:Np);
    for k=1:Np
        njk (:,:,k)= repmat (asd(k),Nc, Nd);
    end
    lbest = njk.*X + (1-njk).*lbest;
    Jold =J;
    % update Clusters centroids
    V = c0*V + c1*rand*(lbest-X) + c2*rand*(repmat(actualbest,1,1,Np)-X);
    X = X+ V;
end


%% Figure 

scatter(O(:,1),O(:,2));
hold on
scatter(actualbest(:,1),actualbest(:,2),80,'r','filled')
scatter(C(:,1),C(:,2),'g','filled')
legend('Data Vector','PSO g_{best}','K-means')
%scatter(X(:,1,1),X(:,2,1),'r*')
%scatter(X0(:,1,1),X0(:,2,1),'k+')