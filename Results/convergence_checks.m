clc
clear all

clusters=2
covariates_habits=3
habits=6
types=3
educ=3
covariates=8
variables_tr=clusters^2*covariates*types
variables_gma=covariates_habits*habits*types
variables_LE=types*2*(clusters+1)*educ;
variables_p=12
generations=37
initial_age=26

cd('C:\Users\jbueren\OneDrive - Istituto Universitario Europeo\endo_health')

fileID=fopen('c_tr.txt');
c_tr=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
c_tr=reshape(c_tr{1},variables_tr,size(c_tr{1},1)/(variables_tr));
iterations=size(c_tr,2)

fileID=fopen('LE.txt');
LE=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
LE=reshape(LE{1},types,2,educ,clusters+1,size(LE{1},1)/(types*2*educ*(clusters+1)));
iterations=size(LE,5)

fileID=fopen('fraction_t.txt');
fraction_t=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
fraction_t=reshape(fraction_t{1},types,2,educ,size(fraction_t{1},1)/(educ*types*2));
iterations=size(fraction_t,4)

fileID=fopen('c_habits.txt');
c_gma=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
c_gma=reshape(c_gma{1},variables_gma,size(c_gma{1},1)/(variables_gma));
iterations=size(c_gma,2)


iterations=min([size(c_gma,2) size(c_tr,2)])

burn=100



%% Histogram from distribution of variables governing transitions

c_l=2
c_l2=1
ind=c_l+(c_l2-1)*c_l2;
figure(2)
for c_l=1:covariates
subplot(2,covariates,c_l)
plot(c_tr(ind+clusters^2*(c_l-1),burn:iterations))
grid on
subplot(2,covariates,covariates+c_l)
hist(c_tr(ind+clusters^2*(c_l-1),burn:iterations))
end



%% Histogram from distribution of variables governing habits

c_l=6 %habits
c_l2=3 %types

ind=1+(c_l-1)*covariates_habits+(c_l2-1)*(covariates_habits*habits);
figure(3)
subplot(2,covariates_habits,1)
plot(c_gma(ind,burn:iterations))
grid on
title('constant')
subplot(2,covariates_habits,2)
plot(c_gma(ind+1,burn:iterations))
grid on
title('age')
subplot(2,covariates_habits,3)
plot(c_gma(ind+2,burn:iterations))
grid on
title('age2')
subplot(2,covariates_habits,4)
hist(c_gma(ind,burn:iterations))
subplot(2,covariates_habits,5)
hist(c_gma(ind+1,burn:iterations))
subplot(2,covariates_habits,6)
hist(c_gma(ind+2,burn:iterations))


%%
max=100
alphas=zeros(habits,generations,types,max);
for it=1:max
    it
gamma=reshape(c_gma(:,iterations-it),covariates_habits,habits,types);
for e_l=1:types
    for h_l=1:habits; for g_l=1:generations
        age=initial_age+(g_l-1)*2-70;
        x(:,1)=[1.0,age,age^2-1];
        alphas(h_l,g_l,e_l,it)=(1-normcdf(0,sum(x(:,1).*gamma(:,h_l,e_l)),1))*100;
    end ;end 
end 
end

colors = {  [0.4660    0.6740    0.1880]    [0.8500    0.3250    0.0980] [0.9290    0.6940    0.1250] [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
pattern = {'none' 'o' 's' '^'};
pattern = {  '--' ':' '-'  '-.' '-'};
    FS=8 %font size
figure(7)
set(7,'position',[150    150    550    250])
ind=0
for h_l=[1 4 5 2 3 6]
    ind=ind+1
    f(h_l)=subplot(2,3,ind)
    for p_l=1:types
        plot(50:2:100,squeeze(mean(alphas(h_l,12:end,p_l,:),4)),'Color',colors{p_l},'linewidth',1.5,'linestyle',pattern{p_l})
        hold on
%         plot(50:2:98,squeeze(prctile(alphas(h_l,:,p_l,:),97.5,4)),'--','Color',colors{p_l})
%         plot(50:2:98,squeeze(prctile(alphas(h_l,:,p_l,:),2.5,4)),'--','Color',colors{p_l})
        ylim([-0.1 1])
        
    % hold off
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
        title('Extreme Obesity','FontWeight','normal','fontsize',FS)
    end 
    yticks([0:25:100])
    xlim([49 100])
    xticks([50:10:100])
%     xlim([26 100])
%     xticks([26:10:100])
    set(gcf,'color','w')
    ylim([-5,105])
    xlabel('Age')
    MS=25 %marker size
set(gca,'FontName','Times New Roman','FontSize',FS);
end

% f(4).Position(1) = 0.25;
% f(5).Position(1) = 0.55;

I=legend('Protective','Harmful','Detrimental','Location','northwest','orientation','horizontal')
legend('boxoff')
I.FontSize=FS
newPosition = [0.45 0.93 0.1 0.1];
    newUnits = 'normalized';
    set(I,'Position', newPosition,'Units', newUnits);
grid off
set(gca,'FontName','Times New Roman','FontSize',FS);

print('C:\Users\jbueren\Dropbox\habits\Slides\v1\figures\health_behaviors','-depsc')





%% Plot Life expectancy for the different groups

%males
for ge_l=1:2
for e_l=1:educ
    [sum(fraction_t(:,ge_l,e_l,:),4)/size(fraction_t,4) sum(LE(:,ge_l,e_l,clusters+1,:),5)/size(LE,5) sum(LE(:,ge_l,e_l,1,:),5)/size(LE,5) sum(LE(:,ge_l,e_l,2,:),5)/size(LE,5)]
end
end

% ind=26
% LE=mean(LE(:,burn:end),2)
% LE_table=[LE(ind:-1:ind-1,end); LE(ind+4:-1:ind+3,end);LE(ind+8:-1:ind+7,end);LE(ind+2:-1:ind+1,end); LE(ind+6:-1:ind+5,end);LE(ind+10:-1:ind+9,end)]
% ind=2
% HLE_table=[LE(ind:-1:ind-1,end); LE(ind+4:-1:ind+3,end);LE(ind+8:-1:ind+7,end);LE(ind+2:-1:ind+1,end); LE(ind+6:-1:ind+5,end);LE(ind+10:-1:ind+9,end)]
% ind=14
% ULE_table=[LE(ind:-1:ind-1,end); LE(ind+4:-1:ind+3,end);LE(ind+8:-1:ind+7,end);LE(ind+2:-1:ind+1,end); LE(ind+6:-1:ind+5,end);LE(ind+10:-1:ind+9,end)]
% ind=2
% fraction_t=mean(fraction_t(:,burn:end),2)
% fraction=[fraction_t(ind:-1:ind-1,end); fraction_t(ind+4:-1:ind+3,end);fraction_t(ind+8:-1:ind+7,end);fraction_t(ind+2:-1:ind+1,end); fraction_t(ind+6:-1:ind+5,end);fraction_t(ind+10:-1:ind+9,end)]
% fraction=fraction.*100
% 
% [fraction LE_table HLE_table ULE_table]
