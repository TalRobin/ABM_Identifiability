clear all
close all

load('Inf_res_1.mat');
colors={'#0072BD','#77AC30', '#A2142F', '#7E2F8E','#D95319'}; %(blue, grean, red, purple, orange)
colors2={[49,196,88],[220,149,27],[26,110,220],[158,93,145]};

No.iter=Dits.num_iter;
No.times=Dits.num_times;
Vrange=Dits.Vrange;
day0=Dits.day0;
tmstep=Dits.tmstep;


for j=1:3 %run over repeats of the inference
    temp=load(['Inf_res_',num2str(j),'.mat']); temp=temp.Dits;
    Gs=temp.Gs;
    Bs=temp.Bs;
    
    %get mean of prior:
    Gmean=zeros(1,No.iter+1);
    Gmean(1,1)=mean(Gs(:,1));
    Bmean=zeros(1,No.iter+1);
    Bmean(1,1)=mean(Bs(:,1));
    
    %get 90% CI of prior:
    G95=zeros(1,No.iter+1);
    G95(1,1)=prctile(Gs(:,1),95);
    B95=zeros(1,No.iter+1);
    B95(1,1)=prctile(Bs(:,1),95);
    G05=zeros(1,No.iter+1);
    G05(1,1)=prctile(Gs(:,1),5);
    B05=zeros(1,No.iter+1);
    B05(1,1)=prctile(Bs(:,1),5);
    
    %get mean and 90% CI of posterior for each itteration:
    for n=1:No.iter
        x=n*(No.times+1)+1;
        x2=(n-1)*(No.times+1)+2:n*(No.times+1);
        Gmean(1,n+1)=mean(Gs(:,x));
        Bmean(1,n+1)=mean(Bs(:,x));
        G95(1,n+1)=prctile(Gs(:,x2),95,"all"); G05(1,n+1)=prctile(Gs(:,x2),5,'all');
        B95(1,n+1)=prctile(Bs(:,x2),95,"all"); B05(1,n+1)=prctile(Bs(:,x2),5,"all");

    end

    %collect to a single structure
    BG(j).Bmean=Bmean;
    BG(j).B95=B95;
    BG(j).B05=B05;
    BG(j).Gmean=Gmean;
    BG(j).G95=G95;
    BG(j).G05=G05;
    clear Bmean B95 B05 Gmean G95 G05 temp Bs Gs

end
clear j



figure(1)
tiledlayout(2,2)

%--------------------------------------------------------------------------
nexttile
hold on
for j=1:3
    fill([0:No.iter, No.iter:-1:0],[BG(j).B95, BG(j).B05(end:-1:1)],colors2{j}./255,FaceAlpha=0.2);
end
for j=1:3
plot(0:No.iter,BG(j).Bmean,'Color',colors2{j}./255,'LineWidth',1);
plot(0:No.iter,BG(j).Bmean,'*','Color',colors2{j}./255);
plot(0:0.2:No.iter,BG(j).Bmean(end),'r.','LineWidth',1.5,'Color',colors2{j}./255)

end
xlabel('Iteration'); ylabel('\beta')
set(gca,'FontSize',17,'FontName','Times New Roman' )
text(-4,0.105,'A.','FontSize',17,'FontName','Times New Roman')
hold off

%--------------------------------------------------------------------------
nexttile

hold on
for j=1:3
    fill([0:No.iter, No.iter:-1:0],[BG(j).G95, BG(j).G05(end:-1:1)],colors2{j}./255,FaceAlpha=0.2);
end
for j=1:3
plot(0:No.iter,BG(j).Gmean,'Color',colors2{j}./255,'LineWidth',1);
plot(0:No.iter,BG(j).Gmean,'*','Color',colors2{j}./255);
plot(0:0.2:No.iter,BG(j).Gmean(end),'r.','LineWidth',1.5,'Color',colors2{j}./255)

end
xlabel('Iteration'); ylabel('\gamma')
set(gca,'FontSize',17,'FontName','Times New Roman' )
%set(gcf,'position',[100,100,600,500])
text(-4,0.105,'B.','FontSize',17,'FontName','Times New Roman')

hold off

%=========================================================================

load("Trajs.mat");
Poss=Trajs.Pos;
P95=prctile(Poss,95,2);
P05=prctile(Poss,5,2);
realTraj=Dits.realTraj;

t=day0+(tmstep:tmstep:tmstep*No.times);
t2=[1:No.times, No.times:-1:1]*tmstep+day0;

%--------------------------------------------------------------------------
nexttile([1 2])

fill(t2,[P95; P05(end:-1:1)],colors2{4}./255,FaceAlpha=0.6);
hold on

plot(t,realTraj,'k','LineWidth',3)
plot(t,Poss(:,1),'Color',colors{1},'LineWidth',1.5)
plot(t,Poss(:,5),'Color',colors{4},'LineWidth',1.5)
set(gca,'FontSize',17,'FontName','Times New Roman' )
xlim([tmstep,tmstep*No.times]+day0)
text(day0-100,42,'C.','FontSize',17,'FontName','Times New Roman')
hold off

set(gcf,'position',[-1500,100,1000,700])
figure(1)
