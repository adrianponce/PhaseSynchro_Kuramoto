%
% Calculate Phase Statistics
% The program extracts the phases of band-passed fMRI signals 
% and calculates: The order parameter, the phase-locking matrix (PLV)
% and the distribution of phase differences
% 
% Input Data: time-series, size NxTxns
%             N: nb. of ROIs
%             T: nb. of time steps
%             ns: nb. of subjects or sessions
%
% Functions needed: surrogates.m and Matlab's Signal Processing Toolbox
%
% A. Ponce-Alvarez, 2015
% see: Ponce-Alvarez et al. (2015) PLoS Comput Biol
% http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004100
% The fMRI data is provided in the Supplemental Dataset of the article.
%--------------------------------------------------------------------------


clear all;

% Load data:
%--------------------------------------------------------------------------
nsubj   = 24; % nb. of subjects
nblocks = 2;  % nb. of blocks

ns = nsubj * nblocks;

T = 300; % nb. of time steps
N = 66;  % nb. of ROIs

Data = zeros(T,N,ns);

n = 0;
for i=1:nblocks
    for j=1:nsubj
        file = ['subj' num2str(j) '_block' num2str(i) '.txt'];
        data = textread(file);
        n = n+1;
        Data(:,:,n)=data;
    end
end


%[T,N,ns]=size(Data); % N rois, T time steps, ns sessions


% Design band-pass fiter
% In principle between 0.04-0.07Hz, see Glerean et al. (2012)
%--------------------------------------------------------------------------
 flp = .04;           % lowpass frequency of filter
 fhi = .07;           % highpass
 npts = T;            % total nb of points
 delt = 2;            % sampling interval

k=2;                  % 2nd order butterworth filter
fnq=1/(2*delt);       % Nyquist frequency
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[b,a]=butter(k,Wn);   % construct the filter

%--------------------------------------------------------------------------

% Initialization:
%--------------------------------------------
Time = (5:T-5)*delt;
Tnew = length(5:T-5); % T after removing borders


KOP = zeros(Tnew,ns);
KOP_0 = zeros(Tnew,ns);

% PLV matrix:
PLV=zeros(N);

% All pairwise combinations:
Comb=nchoosek(1:N,2);
S=size(Comb,1);

% delta-phi distribution:
thetas=-pi:0.1:pi;
nth=length(thetas);
PrPhi=zeros(1,nth);

%nb. of surrogates (it can be reduced to speed the computing time)
repsurr=100;

% prepare paramaters for calculating the spectral density
dt=2; % time resolution: 1 frame = 2 sec.
Ts = T*dt; %% define time of interval
freq = (0:T/2-1)/Ts; %% find the corresponding frequency in Hz
nfreqs=length(freq);
PowSpect=zeros(nfreqs,N,ns);

for subj=1:ns
    
  display(sprintf('subj.:%g',subj))   
  timeseries=Data(:,:,subj)';
  
  % filter
  X_bdpass=timeseries;
  for seed=1:N
  xbp = filtfilt(b,a,timeseries(seed,:));    
  X_bdpass(seed,:) = xbp;    % zero phase filter the data
    
    % spectral density of the bandpass data:
    pw = abs(fft(xbp)); %% absolute value of the fft
    pw = pw(1:floor(T/2)).^2/(T/2); %% take the power of positve freq. half
    PowSpect(:,seed,subj) = pw;

  end
  x=X_bdpass-repmat(mean(X_bdpass,2),[1 T]);
  
  % get phases
  Phases=zeros(N,T);
  for seed=1:N
  Xanalytic = hilbert(x(seed,:));
  Phases(seed,:) = angle(Xanalytic);
  end
 
 % Eliminate the first and last 5 time steps 
 %(to avoid border effects due to the Hilbert transform)
 Phases=Phases(:,5:T-5);

 % Kuramoto order parameter:
 %-----------------------------------------------------
 y=zeros(Tnew,1);
 for t=1:Tnew
  ku=sum(complex(cos(Phases(:,t)),sin(Phases(:,t))))/N;
  y(t)=abs(ku);
 end
 
 % Surrogates, phase randomization 
 % (destroys correlations but preserves the power spectrum 
 % of the time series):
 %--------------------------------------------------------------------
 
     Ynull=zeros(repsurr,Tnew);
     ysurr=zeros(1,Tnew);
     FCplv_Randij=zeros(N);

    for rep=1:repsurr  

      RandPhases=zeros(N,T);

      xsurr_bdp=timeseries;  
      for seed=1:N
      xsurr = surrogates(timeseries(seed,:)); % phase randomization     
      xbp = filtfilt(b,a,xsurr);    
      xsurr_bdp(seed,:) = xbp;    % zero phase filter the data    
      end
      xsurr=xsurr_bdp-repmat(mean(xsurr_bdp,2),[1 T]);


      for seed=1:N
          Xanalytic = hilbert(xsurr(seed,:));
          RandPhases(seed,:) = angle(Xanalytic);
      end
      RandPhases=RandPhases(:,5:T-5);
      for t=1:Tnew
         ku=sum(complex(cos(RandPhases(:,t)),sin(RandPhases(:,t))))/N;
         ysurr(t)=abs(ku);
      end
         Ynull(rep,:)=ysurr;      

      % PLV surrogates   
      fcplv_rand=zeros(N);
       for s=1:S
          k=Comb(s,1);
          l=Comb(s,2); 
          phi=RandPhases(k,:);
          phj=RandPhases(l,:);
          fcplv_rand(k,l)=sum(exp(1i*(phi-phj)))/Tnew;
          fcplv_rand(l,k)=fcplv_rand(k,l);
       end

     FCplv_Randij = FCplv_Randij + abs(fcplv_rand)/repsurr;

    end

     clear xsurr timeseries

     Msync0 = mean(Ynull);
     SDsync0 = std(Ynull); 

     KOP(:,subj)=y;
     KOP_0(:,subj)=Msync0;
 %-------------------------------------------------------------------------    

 
 % Get PLV and the distributio of phase differences:
 %-----------------------------------------------------
 
 fcplv=zeros(N);
 phdiff=zeros(S,Tnew);
 
 for s=1:S
     i=Comb(s,1);
     j=Comb(s,2);
     phi=Phases(i,:);
     phj=Phases(j,:);
     fcplv(i,j)=sum(exp(1i*(phi-phj)))/Tnew;
     fcplv(j,i)=fcplv(i,j);
     dph=phi-phj;
     phdiff(s,:)=(-1).^(dph>0).*2.*pi.*(abs(dph)>=pi) + dph;

 end
 
 % PLV matrix averaged over subjects (corrected using surrogates)
 PLV = PLV + (abs(fcplv)-FCplv_Randij)/ns;

 
% Phase difference distribution:

c=hist(phdiff(:),thetas);
p=c/sum(c);

PrPhi = PrPhi + p/ns;
 
 
end

P1=PrPhi(1);
P2=PrPhi(end);
PrPhi(end)=(P2+P1)/2;
PrPhi(1)=(P2+P1)/2;

% Intrinsic frequencies:
Power_Areas=mean(PowSpect,3);
[~,index]=max(Power_Areas);
Intrinsic_freqs = freq(index);


figure

axes('position',[.13 .65 .8 .3])
plot(Time,KOP(:,1),'k','linewidth',2)
hold on
plot(Time,KOP_0(:,1),'r')
xlabel('Time (s)','fontsize',9)
ylabel('Order parameter','fontsize',9)
legend('data 1','surr.')

axes('position',[.13 .15 .22 .28])
plot(thetas,PrPhi,'k','linewidth',2)
xlabel('\Delta\phi','fontsize',9)
ylabel('Probability','fontsize',9)
set(gca,'xlim',[-pi pi])
title('Phase differences distribution','fontsize',9)

axes('position',[.41 .15 .22 .28])
imagesc(PLV)
xlabel('ROI','fontsize',11)
ylabel('ROI','fontsize',11)
set(gca,'xtick',[],'ytick',[])
pos=get(gca,'position');
h=colorbar;
set(h,'position',[pos(1)+pos(3)+0.01 pos(2) 0.015 pos(4)])
set(gca,'position',pos,'fontsize',8)
title('PLV matrix','fontsize',11)

axes('position',[.77 .15 .17 .28])
hist(Intrinsic_freqs,0.04:.0025:0.07)
set(gca,'xlim',[0.04 0.07])
xlabel('intrinsic frequency (Hz)','fontsize',9)
ylabel('count','fontsize',9)

save Exp_Intrinsic_freqs_and_PLV_share Intrinsic_freqs PLV


