clc
clear all
close all

clusters=2
covariates_habits=4
habits=6
types=3
educ=3
genders=2
variables_p=12
generations=37
initial_age=26
cohorts=5

covariates=2+(types-1)*2
variables_tr=clusters*covariates*educ*genders
variables_gma=covariates_habits*habits*types
variables_H=(clusters+1)*(clusters+1)*generations*types*genders*educ

cd('C:\Users\jbueren\OneDrive - Istituto Universitario Europeo\endo_health')

fileID=fopen('c_tr.txt');
c_tr=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
c_tr=reshape(c_tr{1},covariates,clusters,genders,educ,size(c_tr{1},1)/(variables_tr));


fileID=fopen('LE.txt');
LE=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
LE=reshape(LE{1},types,genders,educ,clusters+1,size(LE{1},1)/(types*genders*educ*(clusters+1)));


fileID=fopen('fraction_t.txt');
fraction_t=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
fraction_t=reshape(fraction_t{1},generations,genders,educ,types,cohorts,size(fraction_t{1},1)/(generations*educ*types*genders*cohorts));

fileID=fopen('fraction_h.txt');
fraction_h=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
fraction_h=reshape(fraction_h{1},generations,clusters,genders,educ,types,size(fraction_h{1},1)/(generations*educ*types*genders*clusters));


fileID=fopen('c_habits.txt');
c_gma=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
c_gma=reshape(c_gma{1},covariates_habits,habits,types,size(c_gma{1},1)/(variables_gma));

fileID=fopen('H.txt');
H=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
H=reshape(H{1},clusters+1,clusters+1,generations,types,genders,educ,size(H{1},1)/(variables_H));


iterations=min([size(c_gma,4) size(c_tr,5)])

burn=1



%% Histogram from distribution of variables governing transitions

for c_l=1:2
ge_l=1
e_l=1

figure(c_l)
for cov_l=1:covariates
subplot(2,covariates,cov_l)
plot(squeeze(c_tr(cov_l,c_l,ge_l,e_l,burn:iterations)))
grid on
subplot(2,covariates,covariates+cov_l)
hist(squeeze(c_tr(cov_l,c_l,ge_l,e_l,burn:iterations)))


end
end


%% Histogram from distribution of variables governing habits

h_l=5 %habits
for y_l=1:types
figure(y_l)
for cov_l=1:covariates_habits
subplot(2,covariates_habits,cov_l)
plot(squeeze(c_gma(cov_l,h_l,y_l,burn:iterations)))
grid on
subplot(2,covariates_habits,cov_l+covariates_habits)
hist(squeeze(c_gma(cov_l,h_l,y_l,burn:iterations)))
end
end


%%
max=1
alphas=zeros(habits,generations,types,max);
ini=iterations-1
for it=0:max-1
for e_l=1:types
    for h_l=1:habits; for g_l=1:generations
        clear x;
        age=initial_age+(g_l-1)*2-70;
        x(:,1)=[1,age,age^2-1,0];
        alphas(h_l,g_l,e_l,it+1)=(1-normcdf(0,sum(x(:,1).*c_gma(:,h_l,e_l,it+ini+1)),1))*100;
    end ;end 
end 
end

alphas(1,1:12,:,:)=NaN;
alphas(4,1:12,:,:)=NaN;
alphas(5,1:12,:,:)=NaN;

colors = { [0.4660    0.6740    0.1880]  [0.9290    0.6940    0.1250]    [0.8500    0.3250    0.0980]   [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
pattern = {  '-'  ':' '--' '-.' '-'};
lw=[1.7 3.0 1.5 ]
    FS=11 %font size
figure(10)
set(10,'position',[150    150    750    350])
ind=0
for h_l=[1 4 5 2 3 6]
    ind=ind+1
    f(h_l)=subplot(2,3,ind)
    for p_l=1:types 
        h(p_l)=plot(26:2:98,squeeze(mean(alphas(h_l,:,p_l,:),4)),'Color',colors{p_l},'linewidth',lw(p_l),'linestyle',pattern{p_l})
        hold on
    end

    if h_l==1
    title('Cancer test','FontWeight','normal','fontsize',FS)
    elseif h_l==2
        title('Drinking','FontWeight','normal','fontsize',FS)
    elseif h_l==3
         title('Smoking','FontWeight','normal','fontsize',FS)
    elseif h_l==4
        title('Cholesterol test','FontWeight','normal','fontsize',FS)
    elseif h_l==5
        title('Flu shot','FontWeight','normal','fontsize',FS) 
    elseif h_l==6
        title('Obesity','FontWeight','normal','fontsize',FS)
    end 
    yticks([0:25:100])
    xlim([25 100])
    xticks([25:10:100])
    set(gcf,'color','w')
    ylim([-5,105])
    xlabel('Age')
    MS=25 %marker size
set(gca,'FontName','Times New Roman','FontSize',FS);
end

I=legend('Protective','Detrimental','Harmful','Location','northwest','orientation','horizontal')
legend('boxoff')
I.FontSize=FS
newPosition = [0.45 0.93 0.1 0.1];
    newUnits = 'normalized';
    set(I,'Position', newPosition,'Units', newUnits);
grid off
set(gca,'FontName','Times New Roman','FontSize',FS);

print('C:\Users\jbueren\Dropbox\habits\Slides\v1\figures\health_behaviors','-depsc')
print('C:\Users\jbueren\Dropbox\habits\Draft\metric_model\figures\health_behaviors','-depsc')



%% Plot Life expectancy for the different groups

for ge_l=1:genders
for e_l=1:educ
    if ge_l==1 && e_l==1
        table=[squeeze(mean(fraction_t(12,ge_l,e_l,:,3,burn:end),6)) mean(LE(:,ge_l,e_l,clusters+1,burn:end),5) mean(LE(:,ge_l,e_l,1,burn:end),5) mean(LE(:,ge_l,e_l,2,burn:end),5)]
    else
        table=vertcat(table,...
               [squeeze(mean(fraction_t(12,ge_l,e_l,:,3,burn:end),6)) mean(LE(:,ge_l,e_l,clusters+1,burn:end),5) mean(LE(:,ge_l,e_l,1,burn:end),5) mean(LE(:,ge_l,e_l,2,burn:end),5)])
    end
    
end
end




%% Plot weights
ge_l=1
e_l=1
FS=10
colors = { [0.4660    0.6740    0.1880]  [0.9290    0.6940    0.1250]    [0.8500    0.3250    0.0980]   [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
pattern = {  '-'  '--' ':' '-.' '-'};
lw=[1.7 1.5 2.0]


%select gender
marker= {'o','s','d' }
ge_l=1
figure(6)
set(6,'position',[150    150    750    750])
for e_l=1:3
    subplot(2,2,e_l)
    for p_l=1:types 
        h(p_l)=errorbar(1910:20:1990,squeeze(mean(fraction_t(12,ge_l,e_l,p_l,:,burn:end),6)),2.*squeeze(std(fraction_t(12,ge_l,e_l,p_l,:,burn:end),0,6)),...
            marker{p_l},'MarkerSize',6,'MarkerFaceColor',colors{p_l})
        h(p_l).Color = colors{p_l}
        hold on
        ylim([0 1])
    end    
    xticks([1910:20:1990])
    xlim([1905 1995])
    if e_l==1
        title('dropout')
    elseif e_l==2
        title('highschool')
    elseif e_l==3
        title('college')
    end
    set(gca,'FontName','Times New Roman','FontSize',FS);
end
hold off
set(gcf,'color','w')
I=legend('Protective','Detrimental','Harmful','Location','northwest','orientation','horizontal')
legend('boxoff')
I.FontSize=FS+1
newPosition = [0.45 0.001 0.1 0.07];
    newUnits = 'normalized';
    set(I,'Position', newPosition,'Units', newUnits);
grid off
set(gca,'FontName','Times New Roman','FontSize',FS);
subplot(2,2,4)
for p_l=1:types 

    LE_v=mean(LE(:,ge_l,3,clusters+1,burn:end),5)'*squeeze(mean(fraction_t(12,ge_l,3,:,:,burn:end),6))-...
        mean(LE(:,ge_l,1,clusters+1,burn:end),5)'*squeeze(mean(fraction_t(12,ge_l,1,:,:,burn:end),6))
    h(p_l)=scatter(1910:20:1990,LE_v,50,"filled",'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7], 'LineWidth',1.5) %,...
    xticks([1910:20:1990])
    ylim([5 12])
    xlim([1905 1995])
    hold on
end
title('life expectancy gradient')
set(gca,'FontName','Times New Roman','FontSize',FS);
print('C:\Users\jbueren\Dropbox\habits\Slides\v1\figures\share_y_cohorts','-depsc')



%% transition pr & fraction by h
ge_l=1
h_l=1
h_l2=3
t_l=1
FS=10
colors = { [0.4660    0.6740    0.1880]  [0.9290    0.6940    0.1250]    [0.8500    0.3250    0.0980]   [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
pattern = {  '-'  '--' ':' '-.' '-'};
lw=[1.7 1.5 2.0]
marker= {'o','s','d' }

figure(7)
set(7,'position',[150    150    750    425])
for e_l=1:3
    subplot(1,3,e_l)
    for p_l=1:types 
        h(p_l)=errorbar(26:4:92,mean(squeeze(H(h_l,h_l2,1:2:34,p_l,ge_l,e_l,burn:end)),2),2.*std(squeeze(H(h_l,h_l2,1:2:34,p_l,ge_l,e_l,burn:end))'),...
            marker{p_l},'MarkerSize',6,'MarkerFaceColor',colors{p_l})
        h(p_l).Color = colors{p_l}
        hold on 
    end 
    if e_l==1
        title('dropout')
    elseif e_l==2
        title('highschool')
    elseif e_l==3
        title('college')
    end
ylim([0 0.7])
xlim([24 91])
xticks(25:10:95)
set(gca,'FontName','Times New Roman','FontSize',FS);
end
set(gcf,'color','w')
I=legend('Protective','Detrimental','Harmful','Location','northwest','orientation','horizontal')
legend('boxoff')
I.FontSize=FS+1
newPosition = [0.45 0.001 0.1 0.07];
newUnits = 'normalized';
set(I,'Position', newPosition,'Units', newUnits);
grid off

figure(8)
set(8,'position',[150    150    750    425])
for e_l=1:3
    subplot(1,3,e_l)
    for p_l=1:types 
        h(p_l)=errorbar(26:4:92,mean(squeeze(fraction_h(1:2:34,h_l,ge_l,e_l,p_l,burn:end)),2),2.*std(squeeze(fraction_h(1:2:34,h_l,ge_l,e_l,p_l,burn:end))'),...
            marker{p_l},'MarkerSize',6,'MarkerFaceColor',colors{p_l})
        h(p_l).Color = colors{p_l}
        hold on 
    end 
    if e_l==1
        title('dropout')
    elseif e_l==2
        title('highschool')
    elseif e_l==3
        title('college')
    end
ylim([0 1])
xlim([24 91])
xticks(25:10:95)
set(gca,'FontName','Times New Roman','FontSize',FS);
end
set(gcf,'color','w')
I=legend('Protective','Detrimental','Harmful','Location','northwest','orientation','horizontal')
legend('boxoff')
I.FontSize=FS+1
newPosition = [0.45 0.001 0.1 0.07];
newUnits = 'normalized';
set(I,'Position', newPosition,'Units', newUnits);
grid off



%% Moments assets

clear all
close all
cd('C:\Users\jbueren\Google Drive\endo_health\metric_model\Results')


    fileID=fopen('wealth_moments_data.txt');
    mean_wealth=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
    fclose(fileID);

    mean_wealth=reshape(mean_wealth{1},8,37,3,3);

for p=1:3
    

colors = {  [0.4660    0.6740    0.1880]   [0.9290    0.6940    0.1250]    [0.8500    0.3250    0.0980] [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
pattern = { '--'  '-'  ':' '-.' '-'};
lw=[1.5 1.5 2.5]
FS=11
figure(2)
set(2,'position',[950    150    750    750])
for e_l=1:3
subplot(3,3,(p-1)*3+e_l)
for y_l=1:3
%     h(y_l)=plot(50:2:98,mean_wealth(3+p,13:37,y_l,e_l)./1000,'Color',colors{y_l},'linewidth',lw(y_l),'linestyle',pattern{y_l})
h(y_l)=plot(26:2:80,mean_wealth(3+p,1:28,y_l,e_l)./1000,'Color',colors{y_l},'linewidth',lw(y_l),'linestyle',pattern{y_l})
hold on
set(gca,'FontName','Times New Roman','FontSize',FS);
end
if e_l==1
    title('dropout')
%     ylim([0 150])
elseif e_l==2
    title('highschool')
%     ylim([0 600])
else
    title('college')
%     ylim([0 1000])
end
ylim([0 2000])
if e_l==1
    if p==1
        ylabel('P25')
    elseif p==2
        ylabel('P50')
    elseif p==3
        ylabel('P75')
    elseif p==4
        ylabel('Mean')
    end 
end
        
end
if p==1
I=legend('Protective','Detrimental','Harmful','Location','northwest','orientation','horizontal')
legend('boxoff')
I.FontSize=FS
newPosition = [0.45 -0.02 0.1 0.1];
newUnits = 'normalized';
set(I,'Position', newPosition,'Units', newUnits);
end
grid off
grid off
set(gca,'FontName','Times New Roman','FontSize',FS);
set(gcf,'color','w')
if p==2
print('C:\Users\jbueren\Dropbox\habits\Slides\v1\figures\wealth_moments_p50','-depsc')
print('C:\Users\jbueren\Dropbox\habits\Draft\metric_model\figures\wealth_moments_p50','-depsc')
elseif p==3
print('C:\Users\jbueren\Dropbox\habits\Slides\v1\figures\wealth_moments_p75','-depsc')
print('C:\Users\jbueren\Dropbox\habits\Draft\metric_model\figures\wealth_moments_p75','-depsc')
end
end
print('C:\Users\jbueren\Dropbox\habits\Slides\v1\figures\wealth_moments','-depsc')
print('C:\Users\jbueren\Dropbox\habits\Draft\metric_model\figures\wealth_moments','-depsc')




%% Labor force participation by educ and health status

clear all
close all
cd('C:\Users\jbueren\Google Drive\endo_health\metric_model\Results')

fileID=fopen('labor_force_participation.txt');
participation=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
participation=reshape(participation{1},5,37,2,3,5);

colors = {  [0.4660    0.6740    0.1880]   [0.9290    0.6940    0.1250]    [0.8500    0.3250    0.0980] [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
pattern = { '--'  '-'  ':' '-.' '-'};
lw=[1.5 1.5 2.5 1.5 2.5]
FS=11
figure(2)
set(2,'position',[150    150    750    350])
y_l=1
c_l=3
for e_l=1:3
subplot(1,3,e_l)
for h_l=1:2
h(h_l)=plot(25:2:65,participation(5,1:21,h_l,e_l,c_l),'Color',colors{h_l},'linewidth',lw(h_l),'linestyle',pattern{h_l})
hold on
set(gca,'FontName','Times New Roman','FontSize',FS);
end
if e_l==1
    title('dropout')
elseif e_l==2
    title('highschool')
else
    title('college')
end
ylim([0.4 1])
xlim([22 68])
xticks(25:5:65)
end
I=legend('Good health','Bad health','Location','northwest','orientation','horizontal')
legend('boxoff')
I.FontSize=FS
newPosition = [0.45 -0.02 0.1 0.1];
newUnits = 'normalized';
set(I,'Position', newPosition,'Units', newUnits);
grid off
grid off
set(gca,'FontName','Times New Roman','FontSize',FS);
set(gcf,'color','w')
print('C:\Users\jbueren\Dropbox\habits\Slides\v1\figures\labor_force_h','-depsc')
print('C:\Users\jbueren\Dropbox\habits\Draft\metric_model\figures\labor_force_h','-depsc')

figure(3)
set(3,'position',[150    150    750    350])
y_l=1
c_l=3
h_l=1
for e_l=1:3
subplot(1,3,e_l)
for c_l=3:5
h(c_l)=plot(25:2:65,participation(5,1:21,h_l,e_l,c_l),'Color',colors{c_l},'linewidth',lw(c_l),'linestyle',pattern{c_l})
hold on
set(gca,'FontName','Times New Roman','FontSize',FS);
end
if e_l==1
    title('dropout')
elseif e_l==2
    title('highschool')
else
    title('college')
end
ylim([0.4 1])
xlim([22 68])
xticks(25:5:65)
end
I=legend('cohort 1950','cohort 1970','cohort 1990','Location','northwest','orientation','horizontal')
legend('boxoff')
I.FontSize=FS
newPosition = [0.45 -0.02 0.1 0.1];
newUnits = 'normalized';
set(I,'Position', newPosition,'Units', newUnits);
grid off
grid off
set(gca,'FontName','Times New Roman','FontSize',FS);
set(gcf,'color','w')
print('C:\Users\jbueren\Dropbox\habits\Slides\v1\figures\labor_force_c','-depsc')
print('C:\Users\jbueren\Dropbox\habits\Draft\metric_model\figures\labor_force_c','-depsc')

%% Labor force participation by educ and health status & previous labor force

clear all
% close all
cd('C:\Users\jbueren\Google Drive\endo_health\metric_model\Results')

fileID=fopen('labor_force_participation_dynamic.txt');
participation=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
participation=reshape(participation{1},5,37,2,3,5,2);

colors = {  [0.4660    0.6740    0.1880]   [0.9290    0.6940    0.1250]    [0.8500    0.3250    0.0980] [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
pattern = { '--'  '-'  ':' '-.' '-'};
lw=[1.5 1.5 2.5 1.5 2.5]
FS=11
h_l=1
for f_l=1:2
figure(f_l+3)
set(f_l+3,'position',[150    150    750    350])
for e_l=1:3
subplot(1,3,e_l)
for c_l=3:5
h(c_l)=plot(26:2:60,participation(5,1:18,h_l,e_l,c_l,f_l),'Color',colors{c_l},'linewidth',lw(c_l),'linestyle',pattern{c_l})
hold on
set(gca,'FontName','Times New Roman','FontSize',FS);
end
if e_l==1
    title('dropout')
elseif e_l==2
    title('highschool')
else
    title('college')
end
ylim([0 1])
xticks(0:10:60)
end
I=legend('cohort 1950','cohort 1970','cohort 1990','Location','northwest','orientation','horizontal')
legend('boxoff')
I.FontSize=FS
newPosition = [0.45 -0.02 0.1 0.1];
newUnits = 'normalized';
set(I,'Position', newPosition,'Units', newUnits);
grid off
grid off
set(gca,'FontName','Times New Roman','FontSize',FS);
set(gcf,'color','w')
end

%% Mean income

clear all
cd('C:\Users\jbueren\Google Drive\endo_health\metric_model\Results')

fileID=fopen('median_income.txt');
mean_wealth=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
mean_wealth=reshape(mean_wealth{1},7,37,3,3,5);

colors = {  [0.4660    0.6740    0.1880]   [0.9290    0.6940    0.1250]    [0.8500    0.3250    0.0980] [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
pattern = { '--'  '-'  ':' '-.' '-'};
lw=[1.5 1.5 2.5]
FS=11
figure(5)
set(5,'position',[150    150    750    350])
for e_l=1:3
subplot(1,3,e_l)
for y_l=1:3
h(y_l)=plot(26:2:64,mean_wealth(7,1:20,y_l,e_l,4)./1000,'Color',colors{y_l},'linewidth',lw(y_l),'linestyle',pattern{y_l})
hold on
set(gca,'FontName','Times New Roman','FontSize',FS);
end
if e_l==1
    title('dropout')
elseif e_l==2
    title('highschool')
else
    title('college')
end
ylim([0 200])
% ylim([10 11.5])
xticks(0:10:60)
% yticks(0:25:150)
end
I=legend('Protective','Detrimental','Harmful','Location','northwest','orientation','horizontal')
legend('boxoff')
I.FontSize=FS
newPosition = [0.45 -0.02 0.1 0.1];
newUnits = 'normalized';
set(I,'Position', newPosition,'Units', newUnits);
grid off
grid off
set(gca,'FontName','Times New Roman','FontSize',FS);
set(gcf,'color','w')
print('C:\Users\jbueren\Dropbox\habits\Slides\v1\figures\mean_income_y','-depsc')
print('C:\Users\jbueren\Dropbox\habits\Draft\metric_model\figures\income_p50','-depsc')

%% Mean across cohorts

clear all
cd('C:\Users\jbueren\Google Drive\endo_health\metric_model\Results')

fileID=fopen('median_income.txt');
mean_wealth=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
mean_wealth=reshape(mean_wealth{1},7,37,3,3,5);

colors = {  [0.4660    0.6740    0.1880]   [0.9290    0.6940    0.1250]    [0.8500    0.3250    0.0980] [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
pattern = { '--'  '-'  ':' '-' '-.' };
lw=[2.5 2.5 2.5 2.0 2.5]
FS=11
figure(6)
set(6,'position',[150    150    750    700])
for p_l=1:2
for e_l=1:3
subplot(2,3,e_l+(p_l-1)*3)
for c_l=[3 4 5]
    stat1=mean_wealth(5,1:20,1,e_l,c_l)./1000
    stat2=mean_wealth(6,1:20,1,e_l,c_l)
    if p_l==1
        h(c_l)=plot(26:2:64,stat1,'Color',colors{c_l},'linewidth',lw(c_l),'linestyle',pattern{c_l})
    else
        h(c_l)=plot(26:2:64,stat2,'Color',colors{c_l},'linewidth',lw(c_l),'linestyle',pattern{c_l})
    end
hold on
if e_l==1 && p_l==1
    ylabel('mean income')
elseif e_l==1 && p_l==2
    ylabel('std log income')
end
set(gca,'FontName','Times New Roman','FontSize',FS);
end
if e_l==1
    title('dropout')
elseif e_l==2
    title('highschool')
else
    title('college')
end
ylim([0 1.5])
% yticks(0:0.2:1.4)
xlim([26 64])
if p_l==1
    ylim([0 180])
%     yticks(0:40:200)
end
end
end

I=legend('cohort 1950','cohort 1970','cohort 1990','Location','northwest','orientation','horizontal')
legend('boxoff')
I.FontSize=FS
newPosition = [0.45 -0.02 0.1 0.1];
newUnits = 'normalized';
set(I,'Position', newPosition,'Units', newUnits);
grid off
grid off
set(gca,'FontName','Times New Roman','FontSize',FS);
set(gcf,'color','w')
print('C:\Users\jbueren\Dropbox\habits\Slides\v1\figures\mean_income_c','-depsc')


