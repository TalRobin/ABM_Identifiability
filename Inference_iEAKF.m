function [Dits] = Inference_iEAKF(Vrange,No,Wantstd,Days,RPos)
%No fields:
% 1. tmstep - number of days in a time step
% 2. Days - number of days in the data set
% 3. vars - number of variables
% 4. Iter - number of itterations
% 5. Ens - number of ensmble members
% 6. IntP - probability to be colonized on day 0

VarNames=Vrange.Properties.VariableNames;

%timesteps: ---------------------------------------------------------------
ts=(1:No.tmstep:No.Days)';
No.times=length(ts)-1;
%--------------------------------------------------------------------------

% Initialization: =========================================================
theta=zeros(No.vars,No.Iter+1);
theta(:,1)=mean(Vrange{:,VarNames},1);
%==========================================================================

% intialize parameters from uniform distribution over the parameter range:=
So=rand(No.ens, No.vars);
So=So.*(ones(No.ens,1)*(Vrange{2,:}-Vrange{1,:}))+(ones(No.ens,1)*Vrange{1,:});
So=array2table(So,'VariableNames',VarNames);

Bs=So{:,'Beta'};
Gs=So{:,'Gamma'};
Poss=So{:,'Beta'};
%==========================================================================

% Iterations: =============================================================
    for n=1:No.Iter
        tic()

        xpost=zeros(No.vars,No.ens,No.times);
        %initialize patient status:
        P_status=rand(No.Pat,No.ens)<No.IntP;
        
        for t=1:No.times
            [P_status,Pos,~,~] = Progress_patients_NEC(Days(ts(t):ts(t+1)),P_status,So);
            
            Poss=[Poss, Pos'];


            % Kalman filter: ==============================================
            Obs=Pos;
            Obs(isnan(Obs))=0;
            Obstruth=RPos(t);

            So=iEAKF(Obs,Obstruth,So,Vrange);

      %--------------------------------------------------------------------
%             Jittering between time steps

            for j=1:2
                if std(So{:,j})<Wantstd(n,j)
                    sig=sqrt((Wantstd(n,j)^2-std(So{:,j})^2));
                    So{:,j}=So{:,j}+normrnd(0,sig,[No.ens,1]);
                    So{So{:,j}<Vrange{1,j},j}=Vrange{1,j}; %check lower bound
                    So{So{:,j}>Vrange{2,j},j}=Vrange{2,j}; %check upper bound
                end
            end

            %--------------------------------------------------------------

            xpost(:,:,t) = (So{:,:})';
            %tscore(:,t)=abs(Obs-obstruth)/obstruth;

            
        end % of time loop
        
        theta(:,n+1)=squeeze(mean(xpost,[2,3]));

        %save the position of each ensmble member at each time point:
        Bs=[Bs, squeeze(xpost(1,:,:))];
        Gs=[Gs, squeeze(xpost(2,:,:))];
        
        %save the posterior ensmble:
        So=squeeze(mean(xpost,3))';
        Bs=[Bs, So(:,1)];
        Gs=[Gs, So(:,2)];
 
        So=array2table(So,'VariableNames',VarNames);
        

        disp(['itteration ', num2str(n),' - ', num2str(toc/60),' mins'])
        
    end % of itteration loop

    Dits.realTraj=RPos;
    Dits.Gs=Gs; 
    Dits.Bs=Bs;
    Dits.Poss=Poss;
    Dits.Vrange=Vrange;
    Dits.tmstep=No.tmstep;
    Dits.IntP=No.IntP;
    Dits.day0=No.day0;
    Dits.num_times=No.times;
    Dits.num_ens=No.ens;
    Dits.num_iter=No.Iter;


end