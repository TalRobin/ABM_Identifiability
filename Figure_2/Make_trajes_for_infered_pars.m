clear all
close all

addpath(fileparts(pwd))

load('Inf_res_1.mat'); %load the ditails of the inference (Dits)
tmstep=Dits.tmstep;
IntP=Dits.IntP;
InBeta=mean(Dits.Bs(:,end)); % beta estiemate from inference
InGamma=mean(Dits.Gs(:,end)); %gamma esteimate frim inference
clear Dits

num_ens=200; %number of trajectories to produce


%load network =============================================================
load('../NY_Network_1.mat');
Days=Network.Days; % Dayly natwork structure
NPat=Network.NPat; %total number of patients
NDays=Network.NDays; %number of days
day0=Network.day0; %date of day before start

RPos=cumsum(Network.daypos.positives);
RPos=diff(RPos(1:tmstep:end));
clear Network
%==========================================================================

VarNames={'Beta','Gamma','Alpha','Rho'};
num_var=size(VarNames,2);
% for a single parameter set: =============================================
Vars=array2table(zeros(1,num_var),'VariableNames',VarNames);
Vars.Beta=InBeta; %Baseline transmission rate, Per day
Vars.Gamma=InGamma; %Importation rate, Per admission
Vars.Alpha(:)=1.5/365; %Patient decolonization rate, Per day
Vars.Rho(:)=0.0160; %Observation rate, Per day
%==========================================================================

%timesteps: ===============================================================
ts=(1:tmstep:NDays)';
num_times=length(ts)-1;
%==========================================================================

P_status=rand(NPat,num_ens)<IntP;
Pos=zeros(num_times,num_ens);

for t=1:num_times
        [P_status,Pos(t,:),~,~] = Progress_patients_uniform(Days(ts(t):ts(t+1)),P_status,Vars);
end


Trajs.Pos=Pos;
Trajs.Vars=Vars;
save("Trajs.mat",'Trajs');
