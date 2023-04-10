function [P_status,Pos,New,Imp] = Progress_patients_NEC(days,P_Pre,Vars)
%Progressing the status of E ensembles by one day
% No enviromental contemination
% returns the positives, the total new patients and the importations

Pos=0;
New=0;
Imp=0;
active_P=days(1).patients;
active_W=days(1).wards;

    for d=2:size(days,2)                                              
                                                                       
        P_status=P_Pre.*(1-ones(size(P_Pre,1),1)*Vars.Alpha');        

        for k=1:size(active_W,1)
                temp=active_P.ward==active_W.ward(k);
                X=sum(P_Pre(active_P.MRN(temp),:).*(active_P.weight(temp)*ones(1,size(P_Pre,2))),1);
                X=X./(active_W.Size(k)); %X=sum(Ck)/(nrk)

                temp2=active_P.MRN(temp);
                temp3=(Vars.Beta').*X; %Beta*X+eps
                temp3=(ones(length(temp2),1)*temp3).*(1-P_Pre(temp2,:));
                P_status(temp2,:)=P_status(temp2,:)+temp3;

        end



        active_P=days(d).patients;
        active_W=days(d).wards;
        new_comer=days(d).Pfirst;
        P_status(new_comer,:)=ones(size(new_comer))*(Vars.Gamma');

        temp=find(~P_status==0);
        P_status(temp)=rand(length(temp),1)<=P_status(temp);

        New=New+sum(P_status>P_Pre);
        Imp=Imp+sum(P_status(new_comer,:));


        P_Pre=P_status;
        
        temp=P_status(active_P.MRN,:);
        test=ones(size(temp,1),1)*(Vars.Rho');
        test=test>rand(size(test));
        Pos=Pos+sum(temp.*test,1);
        %Neg=Neg+sum((1-temp).*test,1);

        clear temp temp2

        
    end
    
    


end