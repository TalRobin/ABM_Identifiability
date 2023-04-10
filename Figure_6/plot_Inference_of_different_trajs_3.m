clear all
close all


load('../Figure_2/Trajs.mat')

Pick=[3,5,1];
Pick=[2,5,1];
%Pick=[2,3,4];

for j=1:length(Pick)
temp=load(['DitsTraj_',num2str(Pick(j)),'b.mat']);
temp.Dits.TTraj=Trajs.Pos(:,j);
Dits(j)=temp.Dits;
end
clear temp j

num_times=Dits(1).num_times;
Vrange=Dits(1).Vrange;
Iter=Dits(1).num_iter;
VarNames={'\beta','\gamma','\alpha','\rho'};
day0=Dits(1).day0;
realTraj=Dits.realTraj;

colors={'#0072BD','#D95319','#77AC30', '#A2142F', '#7E2F8E'}; %(blue, orange, grean, red, purple)
colors2={[26,110,220],[220,149,27],[49,196,88],[158,93,145]};
shapes=['o','d','^','h'];


figure(5)
tiledlayout(2,2)
%==========================================================================
t=(1:num_times)*14+day0;
%figure(52)
nexttile([1,2])
hold on
plot(t,Dits(1).realTraj,'k',LineWidth=2);
for j=1:length(Pick)
    plot(t,Dits(j).TTraj,Color=colors{j},LineWidth=2);
end
xlabel('date'); ylabel('MRSA positives')
set(gca,'FontSize',17,'FontName','Times New Roman' )

x=min(t)-(max(t)-min(t))*0.125; x.Format=t.Format;
text(x,42,'A.','FontSize',17,'FontName','Times New Roman'); 

text(x,-14,'B.','FontSize',17,'FontName','Times New Roman'); 

x=min(t)+(max(t)-min(t))*0.49; x.Format=t.Format;
text(x,-14,'C.','FontSize',17,'FontName','Times New Roman'); 
%==========================================================================

for j=1:width(Dits)
    temp=Dits(j).Bs(:,1:num_times+1:end);
    Dits(j).theta=mean(temp,1);
    Dits(j).theta95=prctile(temp,95);
    Dits(j).theta05=prctile(temp,5);
    
    temp=Dits(j).Gs(:,1:num_times+1:end);
    Dits(j).theta=[Dits(j).theta; mean(temp,1)];
    Dits(j).theta95=[Dits(j).theta95; prctile(temp,95)];
    Dits(j).theta05=[Dits(j).theta05; prctile(temp,5)];

end

%==========================================================================
nexttile
hold on
for j=1:3
y=Dits(j).Bs(:,end); x=Dits(j).Gs(:,end);
fig=scatter(x,y,'MarkerFaceColor',colors{j},'MarkerEdgeColor','k');
fig.MarkerFaceAlpha=0.3;
fig.MarkerEdgeAlpha=0.1;
%fig=scatter(mean(x),mean(y),'MarkerFaceColor',colors{j},'MarkerEdgeColor','k',Marker='x');
%plot(mean(x),mean(y),'Ok','MarkerSize',10,'LineWidth',2)

end
for j=1:3
    y=Dits(j).Bs(:,end); x=Dits(j).Gs(:,end);
    plot(mean(x),mean(y),[shapes(j),'k'],'MarkerSize',8,'LineWidth',2,MarkerFaceColor=colors{j})
end
plot(Trajs.Vars.Gamma,Trajs.Vars.Beta,[shapes(4),'k'],'MarkerSize',15,'LineWidth',2)
Vrange=Dits(1).Vrange;
% xlim([Vrange.Beta(1),Vrange.Beta(2)])
% ylim([Vrange.Gamma(1),Vrange.Gamma(2)])
xlabel('\gamma');ylabel('\beta');
legend('','','','T_1','T_2','T_3','data')
set(gca,'FontSize',17,'FontName','Times New Roman' )


%==========================================================================
% for k=1:2
% nexttile
% hold on
% %plot(1:Iter, realVars{1,k}*ones(1,Iter),'k:')
% for j=1:3
% 
%     theta=Dits(j).theta(k,:);
%     Up=Dits(j).theta95(k,:)-theta;
%     Down=theta-Dits(j).theta05(k,:);
%     y=Dits(j).theta05(k,:);
%     y=[Dits(j).theta95(k,:) , y(end:-1:1)];
%     x=[0:Iter, Iter:-1:0];
%     %errorbar((0:Iter)'+(j-1)*0.1,theta,Up,Down,LineWidth=0.7,Color=colors{j})
%     plot((0:Iter)'+(j-1)*0.1,theta,LineWidth=2,Color=colors{j})
%     fill(x,y,colors2{j}./255,FaceAlpha=0.2);
% 
%     plot((0:Iter)'+(j-1)*0.1,theta,LineWidth=2,Color=colors{j})
% 
% end
% plot(0:Iter,Trajs.Vars{:,k}.*ones(1,Iter+1),'k')
% 
% ylim([0,Vrange{2,k}]);
% xlim([1,Iter+0.3]);
% xlabel('Iteration'); ylabel(VarNames{k})
% hold off
% set(gca,'FontSize',17,'FontName','Times New Roman' )
% 
% end
% %==========================================================================

load('Scan_Infs.mat')

for j=1:height(Scan)
    Pos=Scan.Pos{j};
    mu=mean(Pos,2);
    sig=std(Pos')';
    Pos=Scan.synt_traj{j};
    temp=ones(1,size(Pos,2));
    like=(Pos-mu*temp+0.5)./(sig*temp)./sqrt(2);
    like=-0.5*(erfc(like)-erfc(like-1./(sig*temp)./sqrt(2)));
    like=sum(log(like));

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


LikeTable=zeros(10);
for j=1:10
    for k=1:10
    Pos=Scan.Pos{j};
    mu=mean(Pos,2);
    sig=std(Pos')';
    straj=Scan.synt_traj{k};
    like=(straj-mu+0.5)./(sig)./sqrt(2);
    like=-0.5*(erfc(like)-erfc(like-1./(sig)./sqrt(2)));
    like=sum(log(like));
    LikeTable(j,k)=like;

    like=(realTraj-mu+0.5)./(sig)./sqrt(2);
    like=-0.5*(erfc(like)-erfc(like-1./(sig)./sqrt(2)));
    like=sum(log(like));
    LikeReal(j)=like;
    
    end
end


%Scn.Gamma
nexttile 
hold on
for j=1:length(Pick)
    Ind=Pick(j);
    Scn=Scan(Ind,:);
    errorbar(j,Scn.like_med,Scn.like_med-Scn.like05,Scn.like95-Scn.like_med...
        ,LineWidth=3,Color=colors{j});
    %plot(Scn.Gamma,LikeTable(Pick(j),Pick(j)),'rO',MarkerSize=7,LineWidth=3,Color=colors{j})
    plot(j,LikeReal(Ind),shapes(4),MarkerSize=7,Color='k',LineWidth=3);

    for k=1:length(Pick)

    plot(j,LikeTable(Pick(j),Pick(k)),shapes(k),MarkerSize=7,Color=colors{k})
    plot(j,LikeTable(Pick(j),Pick(k)),shapes(k),MarkerSize=7,Color=colors{k},LineWidth=3)
    end
end
xlabel('\theta_j'); ylabel('log(P(T_i|\theta_j))');
set(gca,'FontSize',17,'FontName','Times New Roman' )
legend('','data','','T_1','','T_2','','T_3')
hold off

Gmax=max(Scan(Pick,:).Gamma); Gmin=min(Scan(Pick,:).Gamma);
xlim([0.5,4]); xticks([1,2,3])
set(gcf,'position',[-1500*0,100*0,1200,700])
