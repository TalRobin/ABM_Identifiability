clear all

num_ens=300;
IntP=0.05;
tmstep=14;

Gammas=linspace(0.06,0.03,15);
Betas=linspace(0,0.03,15);



[Gammas,Betas]=meshgrid(Gammas,Betas);
Gammas=Gammas(:); Betas=Betas(:);

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
%load("Scan_line_3.mat")
for v=height(Scan)+1:num_space
    Vars=Vars_bank(v,:);
    if abs(Vars.Beta+Vars.Gamma-0.058)>0.005
        num_ens=100;
    else 
        num_ens=300;
    end
    
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

    save('Scan_narrow_3.mat','Scan')
    disp(['scan ',num2str(v),' of ', num2str(num_space)])

end
