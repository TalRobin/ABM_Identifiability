function [So] = iEAKF(Obs,Obstruth,So,Vrange)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    VarNames=So.Properties.VariableNames;
    VarNum=size(VarNames,2);
    
    obs_var=1+(0.2*Obstruth).^2;
    prior_var = var(Obs);
    post_var = prior_var*obs_var/(prior_var+obs_var);
    prior_var=max(prior_var, 1e-3); %making sure prior_var>0
    prior_mean = mean(Obs);
    post_mean = post_var*(prior_mean/prior_var + Obstruth/obs_var);
    alpha = (obs_var/(obs_var+prior_var)).^0.5;
    dy = post_mean + alpha*((Obs)-prior_mean)-Obs;
    rr=zeros(1,VarNum);
    for j=1:VarNum
        A=cov(So{:,VarNames(j)},Obs);
        rr(j)=A(2,1)/prior_var;
    end
    dx=dy'*rr;
    
    for j=1:VarNum
    temp=VarNames{j};
    So{:,temp}=(So{:,temp}+dx(:,j));
    So{So{:,temp}<Vrange{1,temp},temp}=Vrange{1,temp}; %check lower bound
    So{So{:,temp}>Vrange{2,temp},temp}=Vrange{2,temp}; %check upper bound
    end

end