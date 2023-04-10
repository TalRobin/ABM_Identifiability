function [P_status,Pos,Nos,Imp] = Progress_patients_uniform(days,P_Pre,Vars)
%Progressing the status of an ensbmble where all Vars are the same

Pos=0;
Nos=0;
Imp=0;
active_P=days(1).patients;
active_W=days(1).wards;
% Npat=size(P_Pre,1);
% Nens=size(P_Pre,2);


    for d=2:size(days,2)                                              
        P_status=P_Pre.*(1-Vars.Alpha);


%         x=zeros(Npat,max(active_P.ward));
%         temp=Npat.*(active_P.ward-1)+active_P.MRN;
%         x(temp)=active_P.weight; %xij= weight of patient i in ward j
% 
%         x=(P_Pre'*x)'; %= xij= weighted sum of all patients in ward i at ensmble member j
%         temp=zeros(size(x));
%         temp(active_W.ward,:)=(active_W.Size)*ones(1,Nens);
%         x=x./temp; %= x/ward size-1;
%         x(temp==0)=0;
%         X=zeros(size(P_status)); X(active_P.MRN,:)=x(active_P.ward,:); % for each patient x of the corresponding ward
%         P_status=P_status+X.*(1-P_Pre).*Vars.Beta;
        
        for k=1:size(active_W,1)
                temp=active_P.ward==active_W.ward(k);
                x=sum(P_Pre(active_P.MRN(temp),:).*(active_P.weight(temp)*ones(1,size(P_Pre,2))),1);
                x=x./(active_W.Size(k)); %X=sum(Ck)/(nrk)

                temp2=active_P.MRN(temp);
                temp3=(Vars.Beta').*x; %Beta*X+eps
                temp3=(ones(length(temp2),1)*temp3).*(1-P_Pre(temp2,:));
                P_status(temp2,:)=P_status(temp2,:)+temp3;

        end


        
        new_comer=days(d).Pfirst;
        P_status(new_comer,:)=Vars.Gamma;
        
        temp=find(~P_status==0);
        P_status(temp)=rand(length(temp),1)<=P_status(temp);
        
        temp=active_P.MRN;
        Nos=Nos+sum(P_status(temp,:)>P_Pre(temp,:));
        Imp=Imp+sum(P_status(new_comer,:));

        active_P=days(d).patients;
        active_W=days(d).wards;


        P_Pre=P_status;
        
        temp=P_status(active_P.MRN,:);
        test=ones(size(temp))*(Vars.Rho);
        test=test>rand(size(test));
        Pos=Pos+sum(temp.*test,1);
        %Neg=Neg+sum((1-temp).*test,1);

    end
    
    


end