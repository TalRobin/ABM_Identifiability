clear all
close all

addpath(fileparts(pwd))
load('../Figure_2/Trajs.mat')

No.Iter=10;%number of iterations
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

for Tnum=1:5
    RPos=Trajs.Pos(:,Tnum);
    Dits=Inference_iEAKF(Vrange,No,Wantstd,Days,RPos);
    save(['DitsTraj_',num2str(Tnum),'.mat'],"Dits")
end


% No.Iter=20;%number of iterations
% No.ens=200; %number of ensmble members
% No.tmstep=2*7; %numer of days in a timestep
% No.IntP=0.05; %probability to be colonized upon admition
% 
% 
% %fprintf('load network')%==================================================
% load('NY_Network_1.mat');
% Days=Network.Days;
% No.Pat=Network.NPat; %total number of patients
% No.ward=Network.NWard; %number of participating wards
% No.Days=Network.NDays; %number of days
% No.day0=Network.day0; %date of day before start
% 
% % Traectory of observations from hospitalization data:
% RPos=cumsum(Network.daypos.positives);
% RPos=diff(RPos(1:No.tmstep:end));
% clear Network
% %==========================================================================
% 
% %fprintf('parameters range:')%=============================================
% VarNames={'Beta','Gamma','Alpha','Rho'};
% Vrange = array2table(zeros(2,length(VarNames)),...
%     'VariableNames',VarNames);
% %[lower,upper]
% Vrange.Beta=[0 ; 0.06]; %Baseline transmission rate, Per day
% Vrange.Gamma=[0.0001 ; 0.06]; %Importation rate, Per admission
% Vrange.Alpha=[1.5 ; 1.5]./365; %Patient decolonization rate, Per day
% Vrange.Rho=[0.0160 ; 0.0160]; %Observation rate, Per day
% 
% No.vars=size(VarNames,2);
% %==========================================================================
% 
% %timesteps: ===============================================================
% ts=(1:tmstep:No.Days)';
% No.times=length(ts)-1;
% %==========================================================================
% 
% % Set desired std in each itteration to control jittering:-----------------
% Wantstd=1./linspace(2,20,No.Iter);
% Wantstd=Wantstd'*diff(Vrange{:,:});
% %--------------------------------------------------------------------------
% 
% for Tnum=5:10
% 
% TPos=Trajs.Pos(:,Tnum);
% 
% % Initialization: =========================================================
% 
% theta=zeros(No.vars,No.Iter+1);
% theta(:,1)=mean(Vrange{:,VarNames},1);
% 
% % intialize parameters from uniform distribution over the parameter range:=
% So=rand(No.ens, No.vars);
% So=So.*(ones(No.ens,1)*(Vrange{2,:}-Vrange{1,:}))+(ones(No.ens,1)*Vrange{1,:});
% So=array2table(So,'VariableNames',VarNames);
% 
% Bs=So{:,'Beta'};
% Gs=So{:,'Gamma'};
% Poss=So{:,'Beta'};
% %==========================================================================
% 
% 
% % Iterations: =============================================================
% %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%     for n=1:Iter
%         tic()
%         
%         xpost=zeros(Novar,No.ens,num_times);
%         %tscore=zeros(num_ens,num_times);
%        
%         P_status=rand(No.Pat,No.ens)<IntP;
%         
%         for t=1:num_times
%             [P_status,Pos,~,~] = Progress_patients_NEC(Days(ts(t):ts(t+1)),P_status,So);
%             
%             Poss=[Poss, Pos'];
% 
% 
%             %Kalman Filter: -------------------------------------------
%             Obs=Pos;
%             Obs(isnan(Obs))=0;
%             obstruth=TPos(t);
% 
%             So=iEAKF(Obs,Obstruth,So,Vrange);
% 
%             %--------------------------------------------------------------------
%             %Jittering between time steps
% 
%             for j=1:2
%                 if std(So{:,j})<Wantstd(n,j)
%                     sig=sqrt((Wantstd(n,j)^2-std(So{:,j})^2));
%                     So{:,j}=So{:,j}+normrnd(0,sig,[No.ens,1]);
%                     So{So{:,j}<Vrange{1,j},j}=Vrange{1,j}; %check lower bound
%                     So{So{:,j}>Vrange{2,j},j}=Vrange{2,j}; %check upper bound
%                 end
%             end
%             
%             %--------------------------------------------------------------
% 
%             xpost(:,:,t) = (So{:,:})';
%             
%         end % of time loop
%         
%         theta(:,n+1)=squeeze(mean(xpost,[2,3]));
% 
%         Bs=[Bs, squeeze(xpost(1,:,:))];
%         Gs=[Gs, squeeze(xpost(2,:,:))];
%         
%         
%         So=squeeze(mean(xpost,3))';
%         Bs=[Bs, So(:,1)];
%         Gs=[Gs, So(:,2)];
%  
%         So=array2table(So,'VariableNames',VarNames);
%         
% 
%         disp(['itteration ', num2str(n),' - ', num2str(toc/60),' mins'])
%         
%     end % of itteration loop
% %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% 
% 
% 
% Dits.realTraj=RPos;
% Dits.Gs=Gs; 
% Dits.Bs=Bs;
% Dits.Poss=Poss;
% Dits.Vrange=Vrange;
% Dits.tmstep=tmstep;
% Dits.IntP=IntP;
% Dits.day0=No.day0;
% Dits.num_times=num_times;
% Dits.num_ens=No.ens;
% Dits.num_iter=Iter;
% 
% save(['DitsTraj_',num2str(Tnum),'b.mat'],"Dits")
% end
