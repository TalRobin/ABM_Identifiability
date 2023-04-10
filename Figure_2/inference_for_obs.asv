clear all
close all

addpath(fileparts(pwd))

No.Iter=20;%number of iterations
No.ens=200; %number of ensmble members
No.tmstep=2*7; %numer of days in a timestep
No.IntP=0.05; %probability to be colonized upon admition



%fprintf('load network')%==================================================
load('NY_Network_1.mat');
Days=Network.Days;
No.Pat=Network.NPat; %total number of patients
No.ward=Network.NWard; %number of participating wards
No.Days=Network.NDays; %number of days
No.day0=Network.day0; %date of day before start

% Traectory of observations from hospitalization data:
RPos=cumsum(Network.daypos.positives);
RPos=diff(RPos(1:No.tmstep:end));
clear Network
%==========================================================================


% parameters range: =======================================================
VarNames={'Beta','Gamma','Alpha','Rho'};
Vrange = array2table(zeros(2,length(VarNames)),...
    'VariableNames',VarNames);
%[lower,upper]
Vrange.Beta=[0 ; 0.1]; %Baseline transmission rate, Per day
Vrange.Gamma=[0.0001 ; 0.1]; %Importation rate, Per admission
Vrange.Alpha=[1.5 ; 1.5]./365; %Patient decolonization rate, Per day
Vrange.Rho=[0.0160 ; 0.0160]; %Observation rate, Per day

%Vrange.Alpha=[realVars.Alpha ; realVars.Alpha];
No.vars=size(VarNames,2);
%==========================================================================

% Set desired std in each itteration to control jittering:-----------------
Wantstd=1./linspace(2,20,No.Iter);
Wantstd=Wantstd'*diff(Vrange{:,:});
%--------------------------------------------------------------------------

for Tnum=1:3
    Dits=Inference_iEAKF(Vrange,No,Wantstd,Days,RPos);
    save(['Inf_res_',num2str(Tnum),'.mat'],"Dits")
end


