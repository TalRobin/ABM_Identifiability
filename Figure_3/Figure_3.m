%Plotting figure 3

clear all
close all

%Load datasets: ===========================================================
%Load likelihood from narrow scan 
Tnarrow=load("Scan_narrow_3.mat"); Tnarrow=Tnarrow.Scan;

%Load likelihood from wide scan
Twide=load("Scan_narrow.mat"); Twide=Twide.Scan;

%Load likelihood allong the line
Tline=load("Scan_line_4.mat"); Tline=Tline.Scan;

% unite likelihood datasets
Tunite=[Twide; Tnarrow; Tline];
%==========================================================================


%extrapulate likelihood landscape from all files ==========================
bu=linspace(min(Tunite.Beta),max(Tunite.Beta),50);
gu=linspace(min(Tunite.Gamma),max(Tunite.Gamma),50);
[gu, bu]=meshgrid(gu, bu);
Lu=griddata(Tunite.Gamma,Tunite.Beta,Tunite.like,gu,bu,'natural');
Lu(isnan(Lu))=min(Lu(:));

%normalize the likelihood
NormConst=log(sum(exp(Lu(:))));
Lu=Lu-NormConst;

%remove very unlikely entries to inprove visibility
Lu(Lu<-500)=-500;
%==========================================================================

% wheited line fitting of ridge line ======================================
W=exp(Tunite.like-max(Tunite.like));
W=W./sum(W);
fx=sum(W.*Tunite.Gamma);
fy=sum(W.*Tunite.Beta);
fxy=sum(W.*Tunite.Beta.*Tunite.Gamma);
fxx=sum(W.*Tunite.Gamma.*Tunite.Gamma);
a=(fxy-fx*fy)/(fxx-fx*fx);
b=-(fxy*fx-fxx*fy)/(fxx-fx*fx);
clear W fx fy fxx fxy
%==========================================================================

%==========================================================================
%Plot main figure =========================================================
figure(3)
tiledlayout(1,2) 

%panel A ------------------------------------------------------------------
nexttile
hold on

g2=gu(1,:); 
b2=bu(:,1);
imagesc('XData',g2,'YData',b2,'CData',Lu)
xlabel('\gamma',FontSize=16,FontWeight='bold'); 
ylabel('\beta',FontSize=16,FontWeight='bold');
colormap('turbo')
colorbar

xlim([min(g2),max(g2)]);
ylim([0,max(b2)]);
title('log-likelihood over parameter landscape')

xl=0:0.002:0.06;
plot(xl,b+a.*xl,'w:',LineWidth=2)
text(0.04,0.05,['\beta=',num2str(b,2),num2str(a,3),'\gamma'],...
    'Color','w','FontSize',17,'Rotation',0,'FontName','Times New Roman' )
set(gca,'FontSize',17,'FontName','Times New Roman' )
hold off

%panel b ------------------------------------------------------------------
%Calculate the likelihood distribution for each parameter set
for j=1:height(Tline)
    Pos=Tline.Pos{j};
    mu=mean(Pos,2);
    sig=std(Pos')';
    temp=ones(1,size(Pos,2));
    like=(Pos-mu*temp+0.5)./(sig*temp)./sqrt(2);
    like=-0.5*(erfc(like)-erfc(like-1./(sig*temp)./sqrt(2)));
    like=sum(log(like));
    Tline.like05(j)=prctile(like,5);
    Tline.like95(j)=prctile(like,95);
    Tline.like_med(j)=prctile(like,50);   
    Tline.like_prntile(j)=sum(like<Tline.like(j))/length(like);
end

nexttile
hold on
errorbar(Tline.Gamma,Tline.like_med-NormConst,...
    Tline.like_med-Tline.like05,Tline.like95-Tline.like_med,'r');
plot(Tline.Gamma,Tline.like-NormConst,'b','LineWidth',3)

xlabel(['\gamma (\beta=',num2str(b,2),'-',num2str(-a,3),'\gamma)']); 
ylabel('Log-likelihood');
title('log-likelihood along ridge')
set(gca,'FontSize',17,'FontName','Times New Roman' )
hold off

set(gcf,'position',[0,100,1500,500])

x=-0.01; y=60;
text(x,y,'B.','FontSize',17,'FontName','Times New Roman');

x=-0.09;
text(x,y,'A.','FontSize',17,'FontName','Times New Roman');
%==========================================================================

%Plot inset ===============================================================
figure(32)
hold on

Lu(Lu<-100)=-100;
s=surf(gu(1,:),bu(:,1),Lu);
xlabel('\gamma',FontSize=16,FontWeight='bold'); ylabel('\beta',FontSize=16,FontWeight='bold');
colormap('turbo')
%temp=find(Sl==max(max(Sl)));
%plot(Sg(temp),Sb(temp),'yo')
xlim([0.03,0.06]);
ylim([0,0.03]);
%zlim([max(Lu(:)-150),max(Lu(:))])
xl=0:0.002:max(0.1);
z=max(max(Lu))*ones(1,length(xl));

plot3(xl,b+a.*xl,z,'b-',LineWidth=6)
set(gca,'FontSize',17,'FontName','Times New Roman' )
view(15,3)
hold off

%========================================================================