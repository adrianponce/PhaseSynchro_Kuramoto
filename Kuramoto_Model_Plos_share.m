%
% Kuramoto model
%
% A. Ponce-Alvarez, 2015
% see: Ponce-Alvarez et al. (2015) PLoS Comput Biol
% http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004100
% The fMRI data is provided in the Supplemental Dataset of the article.
%--------------------------------------------------------------------------

clear all

% load structure (DTI):

load Connectome_Matrix C Order

C=C(Order,Order); %order the connectivity matrix
N =size(C,1);    % nb. of nodes
W = C;
Isubdiag=find(triu(ones(N))-eye(N));

% load empirical data:
% (intrinsic frequencies for each node)

load Exp_Intrinsic_freqs_and_PLV_share Intrinsic_freqs PLV

plv_emp=PLV(Isubdiag);


% Fixed parameters of the simulation:
% Everything is in seconds!

T=5000; % sim. length
dt=0.1; % <--- be careful: maybe you should reduce the dt for high-frequency oscillators!
tspan=0:dt:T;
tds=0:20*dt:T;% 0:10*dt:T; %downsampling
Tmax=length(tds);
Tspan=length(tspan);

f_diff=Intrinsic_freqs'; %intrinsic frequencies (estimated from the data)


w = 2*pi*f_diff;

% global coupling:
gs=0.05:.05:2; %.025:0.025:2;
numGs=length(gs);

Ts = (Tmax-1)*dt*10; %% define time of interval

KOPs     = zeros(Tmax-1,numGs); %kuramoto order parameter
meanR    = zeros(1,numGs); %mean kuramoto order parameter
varR     = zeros(1,numGs); % variance of kuramoto order parameter
PLVModel = zeros(N,N,numGs); % Phase-locking matrix (PLV)
PLVfit   = zeros(1,numGs); % fitting PLV model vs PLV empirical
PLVfit95 = zeros(2,numGs); % confidence interval of the fitting

Comb=nchoosek(1:N,2);
Np=size(Comb,1);


kk=0;

for g=gs

display(g)
G=g;

y=zeros(N,Tmax);

phi=2*pi*rand(N,1);

% transient:
% (eliminate the transient dynamics)
pij=zeros(N);
for t=1:5000
    r = repmat(phi,1,N);
    phi_diffs = r'-r;
    phinew = phi + dt*(w + G*sum(W.*sin(phi_diffs),2)) + .2*sqrt(dt)*randn(N,1);
    phi = mod(phinew,2*pi);    
end

% Start storing phases:
tt=0;
for t=1:Tspan
   
    r = repmat(phi,1,N);
    ynew = phi + dt*(w + G*sum(W.*sin(r'-r),2)) +  .2*sqrt(dt)*randn(N,1);
    phi = mod(ynew,2*pi);
    
    if mod(t,20)==0
    tt=tt+1;    
    y(:,tt)=phi;    
    end    
end

% Kuramoto order param.
kop=abs( sum( exp(1i*y) ) )/N;

% Calculate PLV
PLV_mat=zeros(N);
 for t=1:Tmax-1
  for k=1:Np
     i=Comb(k,1);
     j=Comb(k,2);
     phi=y(i,t);
     phj=y(j,t);
     PLV_mat(i,j)=PLV_mat(i,j) + exp(1i*(phi-phj))/(Tmax-1);
     PLV_mat(j,i)=PLV_mat(j,i) + exp(1i*(phi-phj))/(Tmax-1);
  end
 end
 
 PLV_mat=abs(PLV_mat);
 PLV_mat(1:N+1:end)=1;
 
 kk=kk+1;

 KOPs(:,kk)=kop(1:end-1);
 meanR(kk)=mean(kop(1:end-1)); 
 varR(kk)=mean(kop(1:end-1)); 
 PLVModel(:,:,kk)=PLV_mat;
 
 % Fitting PLV matrix:
 
 plv_model=PLV_mat(Isubdiag);
 
 [rc,pv,ci1,ci2]=corrcoef(plv_emp,plv_model);
 PLVfit(kk)=rc(2);
 PLVfit95(1,kk)=ci1(2);
 PLVfit95(2,kk)=ci2(2);

 
end

figure

axes('position',[.13 .65 .25 .25])
plot(gs,meanR,'linewidth',2)
xlabel('glogbal coupling')
ylabel('<R>')
title('mean order parameter')

axes('position',[.5 .65 .25 .25])
plot(gs,PLVfit,'linewidth',2)
hold on
plot(gs,PLVfit95(1,:),'b:')
plot(gs,PLVfit95(2,:),'b:')
xlabel('glogbal coupling')
ylabel('fit')
title('fitting PLV matrix')
set(gca,'ylim',[0 .6])


