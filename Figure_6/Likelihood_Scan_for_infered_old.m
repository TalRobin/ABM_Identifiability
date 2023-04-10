clear all

addpath(fileparts(pwd))

num_ens=200;
IntP=0.05;
tmstep=14;

% Betas=linspace(0.00001,0.0588,50)';
% Gammas=(Betas-0.05884)./(-1.008);


Gammas=zeros(10,1);
Betas=Gammas;
for j=1:length(Gammas)
    load(['DitsTraj_',num2str(j),'b.mat'])
    Gammas(j)=Dits.Gs(end);
    Betas(j)=Dits.Bs(end);
end

num_space=length(Betas);


VarNames={'Beta','Gamma','Alpha','Rho'};
num_var=size(VarNames,2);
% for a single parameter set: =============================================
Vars_bank=array2table(zeros(num_space,num_var),'VariableNames',VarNames);
Vars_bank.Beta=Betas; %Baseline transmission rate, Per day
Vars_bank.Gamma=Gammas; %Importation rate, Per admission
Vars_bank.Alpha(:)=1.5/365; %Patient decolonization rate, Per day
Vars_bank.Rho(:)=0.0160; %Observation rate, Per day
%==========================================================================




%fprintf('load network')%==================================================
load('NY_Network_1.mat');
Days=Network.Days;
NPat=Network.NPat;
wordnum=Network.NWard;
NDays=Network.NDays;
day0=Network.day0; %date can be retrived by day+day0


RPos=cumsum(Network.daypos.positives);
RPos=diff(RPos(1:tmstep:end));
clear Network
%==========================================================================

%timesteps: ===============================================================
ts=(1:tmstep:NDays)';
num_times=length(ts)-1;
%==========================================================================

Scan=table;
%load("Scan_line_2.mat")
for v=height(Scan)+1:num_space
    Vars=Vars_bank(v,:);
    
    P_status=rand(NPat,num_ens)<IntP;
    Pos=zeros(num_times,num_ens);
    
    for t=1:num_times
        tic()
                [P_status,Pos(t,:),~,~] = Progress_patients_uniform(Days(ts(t):ts(t+1)),P_status,Vars);
                disp(['t=',num2str(t), ' duration:', num2str(toc())])
    end
    
    Scan.Beta(v)=Vars.Beta;
    Scan.Gamma(v)=Vars.Gamma;
    Scan.Alpha(v)=Vars.Alpha;
    Scan.Rho(v)=Vars.Rho;
    Scan.Pos(v)={Pos};
    
    mu=mean(Pos,2);
    sig=std(Pos')';
    like=(RPos-mu+0.5)./sig./sqrt(2);
    like=-0.5*(erfc(like)-erfc(like-1./sig./sqrt(2)));
    like(like==0)=min(like(like>0));
    like=sum(log(like));

    Scan.like(v)=like;

    save('Scan_Infs.mat','Scan')
    figure(1)
    hold on
    plot(Scan.Gamma,Scan.like,'b*')
    hold off

end

Scan=sortrows(Scan,"Gamma");

figure(1)
hold on
plot(Scan.Gamma,Scan.like)
hold off

figure(2)
hold on
plot(Scan.Beta,Scan.like)
hold off

for j=1:height(Scan)
    temp=RPos*ones(1,size(Scan.Pos{j},2));
    temp=abs((Scan.Pos{j}-temp)./temp);
    temp=mean(temp,2);
    Scan.fit(j)=mean(temp);
end


for j=1:height(Scan)
    Pos=Scan.Pos{j};
    mu=mean(Pos,2);
    sig=std(Pos')';
    temp=ones(1,size(Pos,2));
    like=(Pos-mu*temp+0.5)./(sig*temp)./sqrt(2);
    like=-0.5*(erfc(like)-erfc(like-1./(sig*temp)./sqrt(2)));
    like=sum(log(like));
    Scan.like05(j)=prctile(like,5);
    Scan.like95(j)=prctile(like,95);
    Scan.like_med(j)=prctile(like,50);

    
    Scan.like_prntile(j)=sum(like<Scan.like(j))/length(like);
end

figure(3)
hold on
errorbar(Scan.Gamma,Scan.like_med,Scan.like_med-Scan.like05,Scan.like95-Scan.like_med,'r');
plot(Scan.Gamma,Scan.like,'b')
xlabel('\gamma'); ylabel('log-likelihood');
title('log-likelihood along ridge')
hold off

figure(4)
hold on
plot(Scan.Gamma,Scan.like_prntile,'b')
xlabel('\gamma'); ylabel('likelihood prentile');
title('log-likelihood along ridge')
hold off
