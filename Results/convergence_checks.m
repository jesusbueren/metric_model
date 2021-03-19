clc
clear all

clusters=2
covariates_habits=3
habits=5
types=2
educ=3
covariates=13
variables_tr=clusters^2*covariates
variables_gma=covariates_habits*habits*types
variables_LE=types*2*(clusters+1)*educ;
variables_p=12
generations=25
initial_age=50

cd('C:\Users\jbueren\OneDrive - Istituto Universitario Europeo\endo_health')

fileID=fopen('c_tr.txt');
c_tr=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
c_tr=reshape(c_tr{1},variables_tr,size(c_tr{1},1)/(variables_tr));
iterations=size(c_tr,2)

fileID=fopen('LE.txt');
LE=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
LE=reshape(LE{1},variables_LE,size(LE{1},1)/(variables_LE));
iterations=size(LE,2)

fileID=fopen('fraction_t.txt');
fraction_t=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
fraction_t=reshape(fraction_t{1},educ*types*2,size(fraction_t{1},1)/(educ*types*2));
iterations=size(fraction_t,2)

fileID=fopen('c_habits.txt');
c_gma=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
c_gma=reshape(c_gma{1},variables_gma,size(c_gma{1},1)/(variables_gma));
iterations=size(c_gma,2)

fileID=fopen('p.txt')
p=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'})
fclose(fileID)
p=reshape(p{1},variables_p,clusters,size(p{1},1)/variables_p/clusters);

posterior=textread('posterior.txt');

iterations=min([size(c_gma,2) size(p,3) size(c_tr,2)])

figure(1)
burn=500000
plot(posterior(burn:end,1))




%% Histogram from distribution of variables governing transitions

c_l=2
c_l2=2
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

c_l=1 %habits
c_l2=2 %types

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
max=500
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

colors = { [0.8500    0.3250    0.0980] [0.4660    0.6740    0.1880] [0.9290    0.6940    0.1250]   };
pattern = {'none' 'o' 's' '^'};
pattern = {'-' '--'};
figure(7)
set(7,'position',[50    150    700    325*0.75*2])
for h_l=1:5
    f(h_l)=subplot(2,3,h_l)
    for p_l=1:types
        plot(50:2:98,squeeze(prctile(alphas(h_l,:,p_l,:),50,4)),'Color',colors{p_l},'linewidth',2,'linestyle',pattern{p_l})
        hold on
%         plot(50:2:98,squeeze(prctile(alphas(h_l,:,p_l,:),97.5,4)),'--','Color',colors{p_l})
%         plot(50:2:98,squeeze(prctile(alphas(h_l,:,p_l,:),2.5,4)),'--','Color',colors{p_l})
        ylim([-0.1 1])
        xlim([48 100])
    % hold off
    end

    if h_l==1
    title('Cancer test')
    elseif h_l==2
        title('Drinking')
    elseif h_l==3
         title('Smoking')
    elseif h_l==4
        title('Cholesterol test')
    elseif h_l==5
        title('Flu shot') 
    end 
    set(gcf,'color','w')
    ylim([-5,105])
    xlabel('Age')
    FS=11 %font size
FS2=10
MS=25 %marker size
set(gca,'FontName','Times New Roman','FontSize',FS2);
end

f(4).Position(1) = 0.25;
f(5).Position(1) = 0.55;

I=legend('Detrimental','Protective','Location','northwest','orientation','horizontal')
legend('boxoff')
I.FontSize=FS
newPosition = [0.45 0.93 0.1 0.1];
    newUnits = 'normalized';
    set(I,'Position', newPosition,'Units', newUnits);
grid off
set(gca,'FontName','Times New Roman','FontSize',FS2);
print('C:\Users\jbueren\Google Drive\endo_health\draft\figures\health_behaviors','-depsc')



%% Histogram from distributions of bernouilli distributions in ADLs

cl=2
v1=1
v2=2
v3=3
v4=4
v5=5
v6=6
figure(4)
subplot(2,6,1)
plot(squeeze(p(v1,cl,burn:iterations))')
grid on
title('Diff walk')
subplot(2,6,2)
plot(squeeze(p(v2,cl,burn:iterations))')
grid on
title('Diff dress')
subplot(2,6,3)
plot(squeeze(p(v3,cl,burn:iterations))')
grid on
title('Diff bath')
subplot(2,6,4)
plot(squeeze(p(v4,cl,burn:iterations))')
grid on
title('Diff eat')
subplot(2,6,5)
plot(squeeze(p(v5,cl,burn:iterations))')
grid on
title('Diff bed')
subplot(2,6,6)
plot(squeeze(p(v6,cl,burn:iterations))')
grid on
title('Diff toilet')
subplot(2,6,7)
hist(squeeze(p(v1,cl,burn:iterations))')
subplot(2,6,8)
hist(squeeze(p(v2,cl,burn:iterations))')
subplot(2,6,9)
hist(squeeze(p(v3,cl,burn:iterations))')
subplot(2,6,10)
hist(squeeze(p(v4,cl,burn:iterations))')
subplot(2,6,11)
hist(squeeze(p(v5,cl,burn:iterations))')
subplot(2,6,12)
hist(squeeze(p(v6,cl,burn:iterations))')

%% Histogram from distributions of bernouilli distributions in iADLs
cl=2
v1=7
v2=8
v3=9
v4=10
v5=11
v6=12

figure(5)
subplot(2,6,1)
plot(squeeze(p(v1,cl,burn:iterations))')
grid on
title('Diff map')
subplot(2,6,2)
plot(squeeze(p(v2,cl,burn:iterations))')
grid on
title('Diff phone')
subplot(2,6,3)
plot(squeeze(p(v3,cl,burn:iterations))')
grid on
title('Diff money')
subplot(2,6,4)
plot(squeeze(p(v4,cl,burn:iterations))')
grid on
title('Diff medic')
subplot(2,6,5)
plot(squeeze(p(v5,cl,burn:iterations))')
grid on
title('Diff shop')
subplot(2,6,6)
plot(squeeze(p(v6,cl,burn:iterations))')
grid on
title('Diff meals')
subplot(2,6,7)
hist(squeeze(p(v1,cl,burn:iterations))')
subplot(2,6,8)
hist(squeeze(p(v2,cl,burn:iterations))')
subplot(2,6,9)
hist(squeeze(p(v3,cl,burn:iterations))')
subplot(2,6,10)
hist(squeeze(p(v4,cl,burn:iterations))')
subplot(2,6,11)
hist(squeeze(p(v5,cl,burn:iterations))')
subplot(2,6,12)
hist(squeeze(p(v6,cl,burn:iterations))')

%% Median of all variables across clusters
FS=11 %font size
FS2=9
MS=25 %marker size
for cl=1:clusters
for vl=1:variables_p
    p50v(cl,vl)=mean(p(vl,cl,burn:iterations))*100;
end
end

colors = {[0    0.4470    0.7410] [0.8500    0.3250    0.0980] [0.9290    0.6940    0.1250] [0.4940    0.1840    0.5560] [0.4660    0.6740    0.1880]};
pattern = {'o' '^' 'x' 's' '>'};

g=figure(6)
set(6,'position',[50    150    700    325*0.75])
for j=1:clusters
    hold on
    scatter(1:variables_p,[p50v(j,2)';p50v(j,6)';p50v(j,3)';p50v(j,5)';p50v(j,1)';p50v(j,4)'; ...
                                            p50v(j,12)';p50v(j,11)';p50v(j,9)';p50v(j,10)';p50v(j,8)';p50v(j,7)';],MS,pattern{j},'MarkerEdgeColor',colors{j},'LineWidth',1.5)
end
hold off

set(gca,'Xtick',1:variables_p,'XTickLabel',{'DRESS','TOILET','BATH',' BED','WALK','EAT','MEALS','SHOP','MONEY','MEDS','PHONE','MAP'})
set(gcf,'color','w')
xlim([0.5,variables_p+0.5])
ylim([-5,105])
box on
grid off
% yrule = ax.YAxis;
% yrule.FontSize = FS;
grid on
I=legend('Good health','Bad health','Location','northwest','orientation','horizontal')
legend('boxoff')
I.FontSize=FS
newPosition = [0.45 0.92 0.1 0.1];
    newUnits = 'normalized';
    set(I,'Position', newPosition,'Units', newUnits);
grid off
set(gca,'FontName','Times New Roman','FontSize',FS2);
print('C:\Users\jbueren\Google Drive\endo_health\draft\figures\health_states','-depsc')

%% Plot Life expectancy for the different groups

ind=26
LE=mean(LE(:,burn:end),2)
LE_table=[LE(ind:-1:ind-1,end); LE(ind+4:-1:ind+3,end);LE(ind+8:-1:ind+7,end);LE(ind+2:-1:ind+1,end); LE(ind+6:-1:ind+5,end);LE(ind+10:-1:ind+9,end)]
ind=2
HLE_table=[LE(ind:-1:ind-1,end); LE(ind+4:-1:ind+3,end);LE(ind+8:-1:ind+7,end);LE(ind+2:-1:ind+1,end); LE(ind+6:-1:ind+5,end);LE(ind+10:-1:ind+9,end)]
ind=14
ULE_table=[LE(ind:-1:ind-1,end); LE(ind+4:-1:ind+3,end);LE(ind+8:-1:ind+7,end);LE(ind+2:-1:ind+1,end); LE(ind+6:-1:ind+5,end);LE(ind+10:-1:ind+9,end)]
ind=2
fraction_t=mean(fraction_t(:,burn:end),2)
fraction=[fraction_t(ind:-1:ind-1,end); fraction_t(ind+4:-1:ind+3,end);fraction_t(ind+8:-1:ind+7,end);fraction_t(ind+2:-1:ind+1,end); fraction_t(ind+6:-1:ind+5,end);fraction_t(ind+10:-1:ind+9,end)]
fraction=fraction.*100

[fraction LE_table HLE_table ULE_table]
