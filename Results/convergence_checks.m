%% Motivivating graphs HRS
clear all
cd('G:\My Drive\endo_health\data')
hrs=readtable('habits_pr_hrs.csv')
psid=readtable('habits_pr_psid.csv')
colors = { [0.4660    0.6740    0.1880]    [0.8500    0.3250    0.0980]  [0.9290    0.6940    0.1250]   [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
pattern = {  '-'  ':' '--' '-.' '-'}
variable_names = {'cancer_test', 'drunk', 'smoken','cholst','flusht','sport_time'}
% New variable names HRS
newVarNames = {'Cancer', 'Drinking', 'Smoking','Cholesterol','Flu Shot','Exercise'};
hrs = renamevars(hrs, variable_names, newVarNames);

% New variable names PSID
variable_names_psid = {'drinksPerDay', 'nCigsDay'}
newVarNames_psid = { 'Drinking', 'Smoking'};
psid = renamevars(psid, variable_names_psid, newVarNames_psid);
psid = renamevars(psid, 'birth', 'rabyear');
FS=10
figure(1)
set(1,'position',[150    150    500    250])
for e_l=[1 3]
    for v_l=1:length(newVarNames)
        var_name = newVarNames{v_l};
        var_name2=strcat("N_",variable_names{v_l})
        
        
        subplot(2,3,v_l)
        for c_l=1:3
            filtered_hrs = hrs(hrs.raeduc == e_l & hrs.rabyear == c_l &  hrs.(var_name2)>50, :)
            plot(filtered_hrs.int_age, filtered_hrs.(var_name).*100,"color",colors{c_l},"LineWidth",2,'linestyle',pattern{e_l})
            if v_l==2 || v_l==3
                var_name2_psid=strcat("N_",variable_names_psid{v_l-1})
                filtered_psid = psid(psid.raeduc == e_l & psid.rabyear == c_l &  psid.(var_name2_psid)>50, :)
                if isempty(filtered_psid.(var_name))==0 
                    plot(filtered_psid.int_age, filtered_psid.(var_name).*100,"color",colors{c_l},"LineWidth",2,'linestyle',pattern{e_l})
                end
            end
            hold on
            title(var_name,'FontWeight','normal')
        end
        
%         ylim([0 100])
    end
end
% legend([p(1,1) p(2,1) p(3,1) p(1,3) p(2,3) p(3,3)],'HSD: 1930','HSD: 1950','HSD: 1970','CG: 1930','CG: 1950','CG: 1970')
        set(gcf,'color','w')
        set(gca,'FontName','Times New Roman','FontSize',FS);

%% Motivivating graphs PSID
clear all
cd('C:\Users\jbueren\Google Drive\endo_health\data')
hrs=readtable('habits_pr_psid.csv')
colors = { [0.4660    0.6740    0.1880]    [0.8500    0.3250    0.0980]  [0.9290    0.6940    0.1250]   [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
pattern = {  '-'  ':' '--' '-.' '-'}
variable_names = {'drinksPerDay', 'nCigsDay'}
% New variable names
newVarNames = { 'Drinking', 'Smoking'};
hrs = renamevars(hrs, variable_names, newVarNames);
hrs = renamevars(hrs, 'birth', 'rabyear');

for e_l=[1 3]
    for v_l=1:length(newVarNames)
        var_name = newVarNames{v_l};
        var_name2=strcat("N_",variable_names{v_l})
        figure(1)
        subplot(1,2,v_l)
        for c_l=1:3
            filtered_hrs = hrs(hrs.raeduc == e_l & hrs.rabyear == c_l &  hrs.(var_name2)>40, :)
            plot(filtered_hrs.int_age, filtered_hrs.(var_name).*100,"color",colors{c_l},"LineWidth",2,'linestyle',pattern{e_l})
            hold on
            title(var_name,'FontWeight','normal')
        end
        ylim([0 100])
    end
end

%% Metric model results
clc
clear all
close all

types=2 % select the number of health behavior groups
clusters=2
covariates_habits=4
habits=3
covariates_habits_med=5
habits_med=3
types_s=num2str(types)
educ=3
genders=2
variables_p=12
generations=37
initial_age=26
cohorts=5

covariates=types*3
variables_tr=clusters*covariates*educ*genders
variables_gma=covariates_habits*habits*types
variables_gma_med=covariates_habits_med*habits_med*types

variables_H=(clusters+1)*(clusters+1)*generations*types*genders*educ
covariates_mixture=cohorts+1
variables_delta=covariates_mixture*(types-1)*genders*educ

cd('C:\Users\Jesus Bueren\results_local\endo_health')

fileID=fopen(strcat('c_tr_d_',types_s,'.txt'));
c_tr=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
c_tr=reshape(c_tr{1},covariates,clusters,genders,educ,size(c_tr{1},1)/(variables_tr));


fileID=fopen(strcat('LE_',types_s,'.txt'));
LE=textscan(fileID,'%20.8f','TreatAsEmpty',{'**************'});
fclose(fileID);
LE=reshape(LE{1},types,genders,educ,clusters+1,size(LE{1},1)/(types*genders*educ*(clusters+1)));


fileID=fopen(strcat('fraction_t_',types_s,'.txt'));
fraction_t=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
fraction_t=reshape(fraction_t{1},generations,genders,educ,types,cohorts,size(fraction_t{1},1)/(generations*educ*types*genders*cohorts));

fileID=fopen(strcat('fraction_h_',types_s,'.txt'));
fraction_h=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
fraction_h=reshape(fraction_h{1},generations,clusters,genders,educ,types,size(fraction_h{1},1)/(generations*educ*types*genders*clusters));


fileID=fopen(strcat('c_habits_',types_s,'.txt'));
c_gma=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
c_gma=reshape(c_gma{1},covariates_habits,habits,types,size(c_gma{1},1)/(variables_gma));

fileID=fopen(strcat('c_habits_med_',types_s,'.txt'));
c_gma_med=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
c_gma_med=reshape(c_gma_med{1},covariates_habits_med,habits,types,size(c_gma_med{1},1)/(variables_gma_med));

fileID=fopen(strcat('H_',types_s,'.txt'));
H=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
H=reshape(H{1},clusters+1,clusters+1,generations,types,genders,educ,size(H{1},1)/(variables_H));

fileID=fopen(strcat('delta_',types_s,'.txt'));
c_delta=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
c_delta=reshape(c_delta{1},covariates_mixture*(types-1),genders,educ,size(c_delta{1},1)/(variables_delta));


iterations=min([size(c_gma,4) size(c_tr,5) size(LE,5)])



burn=1
e_l=2
figure(1)
for t_l=1:types-1
for c_l=1:covariates_mixture
subplot(types-1,covariates_mixture,c_l+(t_l-1)*covariates_mixture )
plot(squeeze(c_delta(c_l,1,e_l,:)))
end
end



figure(2)
for e_l=1:3
subplot(1,3,e_l)
plot(squeeze(fraction_t(:,1,e_l,:,3,end)))
ylim([0 1])
end


%% Histogram from distribution of variables governing transitions

for e_l=1:3
for c_l=1:2
ge_l=1
figure('units','normalized','outerposition',[0 0 1 1])
for cov_l=1:covariates
subplot(2,covariates,cov_l)
plot(squeeze(c_tr(cov_l,c_l,ge_l,e_l,burn:iterations)))
grid on
subplot(2,covariates,covariates+cov_l)
hist(squeeze(c_tr(cov_l,c_l,ge_l,e_l,burn:iterations)))

end
end
end


%% Histogram from distribution of variables governing habits

h_l=1 %habits
for y_l=1:types
figure(y_l)
for cov_l=1:covariates_habits_med
subplot(2,covariates_habits_med,cov_l)
plot(squeeze(c_gma_med(cov_l,h_l,y_l,burn:iterations)))
grid on
subplot(2,covariates_habits_med,cov_l+covariates_habits_med)
hist(squeeze(c_gma_med(cov_l,h_l,y_l,burn:iterations)))
end
end


%%
% close all
max = 1;
alphas = zeros(6, generations, types, max);
ini = iterations - 1;
bh = 0; % bad health dummy
ins=1; % insurance status dummy

ages = initial_age + (0:generations-1) * 2 ;  % Precomputar edades
x = [ones(generations, 1), ages', (ages.^2 - 1)', bh * ones(generations, 1)]; % Matriz de regresores
x_med = [ones(generations, 1), ages', (ages.^2 - 1)', bh * ones(generations, 1), ins * ones(generations, 1)]; % Matriz de regresores

for it = burn:iterations
    it_l = it - burn + 1;
    for e_l = 1:types
        for h_l = 1:3

            c_vec = squeeze(c_gma(:, h_l, e_l, it));  % Extraer coeficientes de una vez
            alphas(h_l, :, e_l, it_l) = (1 - normcdf(0, x * c_vec, 1)) * 100;

            c_vec = squeeze(c_gma_med(:, h_l, e_l, it));  % Extraer coeficientes de una vez
            alphas(h_l+3, :, e_l, it_l) = (1 - normcdf(0, x_med * c_vec, 1)) * 100;

        end
    end
end

ini_v(1:2)=1
ini_v(3:6)=15
alphas(3:6,1:14,:,:)=NaN;

if types==2
    colors = { [0.4660    0.6740    0.1880]    [0.8500    0.3250    0.0980]  [0.9290    0.6940    0.1250]   [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
elseif types==3
    colors = { [0.4660    0.6740    0.1880]     [0.9290    0.6940    0.1250] [0.8500    0.3250    0.0980]   [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
elseif types==4
    colors = { [0.4660    0.6740    0.1880]     [0   0.4470    0.7410]   [0.9290    0.6940    0.1250]  [0.8500    0.3250    0.0980] [0.4940    0.1840    0.5560]};
end

if types==2
    pattern = {  '-'  ':'};
elseif types==3
    pattern = {  '-'   '--' ':'};
elseif types==4
    pattern = {  '-'   '--' '-.' ':'};
end

lw=[1.2 1.5 1.5 1.5 ] %[1.7 2.0 1.5 1.5 ]
    FS=8 %font size
figure(10)
set(10,'position',[150    150    500    250])
ind=0
for h_l=1:6
    ind=ind+1
    f(h_l)=subplot(2,3,ind)
    for p_l=1:types 
        mean_val = squeeze(mean(alphas(h_l,:,p_l,:),4));
        lower_bound = squeeze(prctile(alphas(h_l,:,p_l,:), 2.5, 4));
        upper_bound = squeeze(prctile(alphas(h_l,:,p_l,:), 97.5, 4));
        x_vals = 26:2:98;
        fill([x_vals(ini_v(h_l):end), fliplr(x_vals(ini_v(h_l):end))], [lower_bound(ini_v(h_l):end), fliplr(upper_bound(ini_v(h_l):end))], ...
            colors{p_l}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold on
        h(p_l) = plot(x_vals, mean_val, 'Color', colors{p_l}, ...
            'linewidth', lw(p_l), 'linestyle', pattern{p_l});
        hold on
    end

    if h_l==1
        title('Drinking','FontWeight','normal','fontsize',FS)
    elseif h_l==2
        title('Smoking','FontWeight','normal','fontsize',FS)
    elseif h_l==3
        title('Exercise','FontWeight','normal','fontsize',FS) 
    elseif h_l==4
        title('Cancer test','FontWeight','normal','fontsize',FS)
    elseif h_l==5
        title('Cholesterol test','FontWeight','normal','fontsize',FS)
    elseif h_l==6
        title('Flu shot','FontWeight','normal','fontsize',FS) 
    end 
    yticks(0:25:100)
    xlim([25 100])
    xticks(30:20:100)
    set(gcf,'color','w')
    ylim([-5,105])
    if h_l==1
        yticks(0:5:20)
        ylim([-2,22])
    end
    xlabel('Age')
    MS=25 %marker size
set(gca,'FontName','Times New Roman','FontSize',FS);
end
if types==2
    I=legend([h(1),h(2)],'Protective','Detrimental','Location','northwest','orientation','horizontal')
elseif types==3
    I=legend([h(1),h(2),h(3)],'Group 1','Group 2 ','Group 3','Location','northwest','orientation','horizontal')
elseif types==4
    I=legend([h(1),h(2),h(3),h(4)],'Group 1','Group 2','Group 3','Group 4','Location','northwest','orientation','horizontal')
end 
legend('boxoff')
I.FontSize=FS
newPosition = [0.45 0.93 0.1 0.1];
    newUnits = 'normalized';
    set(I,'Position', newPosition,'Units', newUnits);
grid off
set(gca,'FontName','Times New Roman','FontSize',FS);
if bh==0
    print(strcat('C:\Users\Jesus Bueren\Dropbox\habits\Draft\figures\health_behaviors',types_s),'-depsc')
else
    print(strcat('C:\Users\Jesus Bueren\Dropbox\habits\Draft\figures\health_behaviors',types_s,'_bh'),'-depsc')
end

figure(10)
set(10,'position',[150    150    500    220])
% print('C:\Users\jbueren\Dropbox\habits\Slides\2024_EUI_PhD\figures\health_behaviors','-depsc')



%% Plot Life expectancy for the different groups
close all
if types==2
    colors = { [0.4660    0.6740    0.1880]    [0.8500    0.3250    0.0980]  [0.9290    0.6940    0.1250]   [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
elseif types==3
    colors = { [0.4660    0.6740    0.1880]     [0.9290    0.6940    0.1250] [0.8500    0.3250    0.0980]   [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
elseif types==4
    colors = { [0.4660    0.6740    0.1880]     [0   0.4470    0.7410]   [0.9290    0.6940    0.1250]  [0.8500    0.3250    0.0980] [0.4940    0.1840    0.5560]};
end
lw=[1.2 1.5 1.5 1.5 ] 

figure(1)
for e_l=1:3
    subplot(1,3,e_l)
    for p_l=1:types 
        h(p_l)=plot(squeeze(LE(p_l,1,e_l,clusters+1,burn:iterations))', 'Color', colors{p_l}, ...
            'linewidth', lw(p_l))
        hold on
    end
    if e_l==1
        if types==2
            I=legend([h(1),h(2)],'Protective','Detrimental','Location','northwest','orientation','horizontal')
        elseif types==3
            I=legend([h(1),h(2),h(3)],'Group 1','Group 2', 'Group 3','Location','northwest','orientation','horizontal')
        end
    end
    ylim([15 40])
end
c_l=3
for ge_l=1:1 %genders
for e_l=1:educ
    if ge_l==1 && e_l==1
        table=[ squeeze(mean(fraction_t(12,ge_l,e_l,:,c_l,burn:iterations),6)) mean(LE(:,ge_l,e_l,clusters+1,burn:iterations),5) mean(LE(:,ge_l,e_l,1,burn:end),5) mean(LE(:,ge_l,e_l,2,burn:iterations),5)]
    else
        table=vertcat(table,...
               [ squeeze(mean(fraction_t(12,ge_l,e_l,:,c_l,burn:iterations),6)) mean(LE(:,ge_l,e_l,clusters+1,burn:iterations),5) mean(LE(:,ge_l,e_l,1,burn:iterations),5) mean(LE(:,ge_l,e_l,2,burn:iterations),5)])
    end
    
end
end

for e_l=1:educ
    Av_LE(e_l,:)=sum(squeeze(fraction_t(12,ge_l,e_l,:,c_l,burn:iterations)).*squeeze(LE(:,1,e_l,clusters+1,burn:iterations)),1)
    Av_LE_c(e_l,:)=sum(squeeze(fraction_t(12,ge_l,3,:,c_l,burn:iterations)).*squeeze(LE(:,1,e_l,clusters+1,burn:iterations)),1)
end

figure(2)
subplot(1,3,1)
plot(Av_LE(3,:))
hold on
plot(Av_LE(1,:))
plot(Av_LE_c(1,:))
subplot(1,3,2)
plot(squeeze(fraction_t(12,ge_l,1,1,c_l,burn:iterations)))
hold on
plot(squeeze(fraction_t(12,ge_l,3,1,c_l,burn:iterations)))
subplot(1,3,3)
plot((Av_LE_c(1,:)-Av_LE(1,:))./(Av_LE(3,:)-Av_LE(1,:)))



%% Plot weights across cohorts
e_l=1

if types==2
    colors = { [0.4660    0.6740    0.1880]    [0.8500    0.3250    0.0980]  [0.9290    0.6940    0.1250]   [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
elseif types==3
    colors = { [0.4660    0.6740    0.1880]     [0.9290    0.6940    0.1250] [0.8500    0.3250    0.0980]   [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
elseif types==4
    colors = { [0.4660    0.6740    0.1880]     [0   0.4470    0.7410]   [0.9290    0.6940    0.1250]  [0.8500    0.3250    0.0980] [0.4940    0.1840    0.5560]};
end
pattern = {  '-'  '--' ':' '-.' '-'};
lw=[1.7 1.5 2.0 1.5]

FS=8

%select gender
marker= {'o','s','d' }
ge_l=1
figure(6)
set(6,'position',[150    150    500    280])
for e_l=1:3
    subplot(2,2,e_l)
    for p_l=1:types 
        mean_val2=squeeze(mean(fraction_t(12,ge_l,e_l,p_l,:,burn:iterations), 6));
        lower_bound2 = squeeze(prctile(fraction_t(12,ge_l,e_l,p_l,:,burn:iterations), 2.5, 6));
        upper_bound2 = squeeze(prctile(fraction_t(12,ge_l,e_l,p_l,:,burn:iterations), 97.5, 6));
        x_vals2=10:20:90;
        h(p_l) = plot(x_vals2, mean_val2, 'Color', colors{p_l}, ...
            'linewidth', lw(p_l), 'linestyle', pattern{p_l});
        hold on
        fill([x_vals2, fliplr(x_vals2)], [lower_bound2', fliplr(upper_bound2')], ...
            colors{p_l}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        ylim([0 1])
    end  
    xlabel('Birth Year')
    xticks([10:20:90])
    yticks([0:0.25:1])
    xlim([5 95])
    if e_l==1
        title('HSD','FontWeight','Normal')
    elseif e_l==2
        title('HSG','FontWeight','Normal')
    elseif e_l==3
        title('CG','FontWeight','Normal')
    end
    set(gca,'FontName','Times New Roman','FontSize',FS);
end
hold off
set(gcf,'color','w')
if types==2
    I=legend([h(1),h(2)],'Protective','Detrimental','Location','northwest','orientation','horizontal')
elseif types==3
    I=legend([h(1),h(2),h(3)],'Group 1','Group 2', 'Group 3','Location','northwest','orientation','horizontal')
elseif types==4
    I=legend([h(1),h(2),h(3),h(4)],'Group 1','Group 2', 'Group 3', 'Group 4','Location','northwest','orientation','horizontal')
end

legend('boxoff')
I.FontSize=FS+1
newPosition = [0.45 0.94 0.1 0.07];
    newUnits = 'normalized';
    set(I,'Position', newPosition,'Units', newUnits);
grid off
set(gca,'FontName','Times New Roman','FontSize',FS);
subplot(2,2,4)
for c_l=1:cohorts
LE_v(c_l)=mean(sum(squeeze(fraction_t(12,ge_l,3,:,c_l,burn:iterations)).*squeeze(LE(:,ge_l,3,clusters+1,burn:iterations)))-...
               sum(squeeze(fraction_t(12,ge_l,1,:,c_l,burn:iterations)).*squeeze(LE(:,ge_l,1,clusters+1,burn:iterations))))
lower_bound2(c_l) = squeeze(prctile(sum(squeeze(fraction_t(12,ge_l,3,:,c_l,burn:iterations)).*squeeze(LE(:,ge_l,3,clusters+1,burn:iterations)))-...
               sum(squeeze(fraction_t(12,ge_l,1,:,c_l,burn:iterations)).*squeeze(LE(:,ge_l,1,clusters+1,burn:iterations))), 2.5));
        upper_bound2(c_l) = squeeze(prctile(sum(squeeze(fraction_t(12,ge_l,3,:,c_l,burn:iterations)).*squeeze(LE(:,ge_l,3,clusters+1,burn:iterations)))-...
               sum(squeeze(fraction_t(12,ge_l,1,:,c_l,burn:iterations)).*squeeze(LE(:,ge_l,1,clusters+1,burn:iterations))), 97.5));
end
% h2=scatter(10:20:90,LE_v,10,"filled",'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7], 'LineWidth',1.5) 
h2=plot(x_vals2, LE_v, 'Color', colors{4}, ...
            'linewidth', lw(4));
hold on
fill([x_vals2, fliplr(x_vals2)], [lower_bound2', fliplr(upper_bound2')], ...
            colors{4}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xticks([10:20:90])
% yticks([4:1:9])
%  ylim([4 9])
xlim([05 95])
hold on
xlabel('Birth Year')
title('LE gradient','FontWeight','Normal')
set(gca,'FontName','Times New Roman','FontSize',FS);
% print('C:\Users\jbueren\Dropbox\habits\Draft\figures\share_y_cohorts','-depsc')
set(6,'position',[150    150    500    220])
% print('C:\Users\jbueren\Dropbox\habits\Slides\2024_EUI_PhD\figures\share_y_cohorts','-depsc')
%% Plot weights across age for a given cohort
ge_l=1
e_l=1
c_l=3
FS=10
if types==2
    colors = { [0.4660    0.6740    0.1880]    [0.8500    0.3250    0.0980]  [0.9290    0.6940    0.1250]   [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
elseif types==3
    colors = { [0.4660    0.6740    0.1880]     [0.9290    0.6940    0.1250] [0.8500    0.3250    0.0980]   [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
elseif types==4
    colors = { [0.4660    0.6740    0.1880]     [0   0.4470    0.7410]   [0.9290    0.6940    0.1250]  [0.8500    0.3250    0.0980] [0.4940    0.1840    0.5560]};
end
endpattern = {  '-'  '--' ':' '-.' '-'};
lw=[1.7 1.5 2.0 1.5]


%select gender
marker= {'o','s','d','o' }
ge_l=1
figure(6)
set(6,'position',[150    150    500    250])
% for ge_l=1:2
ge_l=1
for e_l=1:3
    subplot(1,3,e_l)
    for p_l=1:types 
        h(p_l)=errorbar(26:4:100,squeeze(mean(fraction_t(1:2:end,ge_l,e_l,p_l,c_l,burn:iterations),6)).*100,2.*squeeze(std(fraction_t(1:2:end,ge_l,e_l,p_l,c_l,burn:iterations),0,6)).*100,...
            marker{p_l},'MarkerSize',6,'MarkerFaceColor',colors{p_l})
        h(p_l).Color = colors{p_l}
        hold on
        ylim([0 100])
    end    
    yticks([0:25:100])
    xlim([25 100])
    xticks([25:10:100])
    set(gcf,'color','w')
    ylim([-5,105])
    xlabel('Age')
    if e_l==1
        title('dropout','FontWeight','normal','fontsize',FS)
    elseif e_l==2
        title('highschool','FontWeight','normal','fontsize',FS)
    elseif e_l==3
        title('college','FontWeight','normal','fontsize',FS)
    end
    set(gca,'FontName','Times New Roman','FontSize',FS);
end
% end
% hold off
% set(gcf,'color','w')
% I=legend('Protective','Detrimental','Harmful','Location','northwest','orientation','horizontal')
% legend('boxoff')
% I.FontSize=FS+1
% newPosition = [0.45 0.95 0.1 0.07];
%     newUnits = 'normalized';
%     set(I,'Position', newPosition,'Units', newUnits);
% grid off
set(gca,'FontName','Times New Roman','FontSize',FS);

% print('C:\Users\jbueren\Dropbox\habits\Slides\v2\figures\share_y_age','-depsc')



%% transition pr & fraction by h
ge_l=1
FS=11
if types==2
    colors = { [0.4660    0.6740    0.1880]    [0.8500    0.3250    0.0980]  [0.9290    0.6940    0.1250]   [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
elseif types==3
    colors = { [0.4660    0.6740    0.1880]     [0.9290    0.6940    0.1250] [0.8500    0.3250    0.0980]   [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
elseif types==4
    colors = { [0.4660    0.6740    0.1880]     [0   0.4470    0.7410]   [0.9290    0.6940    0.1250]  [0.8500    0.3250    0.0980] [0.4940    0.1840    0.5560]};
end
pattern = {  '-'  '--' ':' '-.' '-'};
lw=[1.7 1.5 2.0]
marker= {'o','s','d','d' }

h_str=["gh","bh"]
h2_str=["gh","bh","D"]

for h_l=1:2
for h_l2=1:3



figure(7)
set(7,'position',[150    100    750    900])
for e_l=1:3
    subplot(6,3,e_l+(h_l2-1)*3+(h_l-1)*9)
    for p_l=1:types 
        h(p_l)=errorbar(26:4:92,mean(squeeze(H(h_l,h_l2,1:2:34,p_l,ge_l,e_l,burn:end)),2),2.*std(squeeze(H(h_l,h_l2,1:2:34,p_l,ge_l,e_l,burn:end))'),...
            marker{p_l},'MarkerSize',6,'MarkerFaceColor',colors{p_l})
%         h(p_l)=errorbar(26:4:92,mean(squeeze(H(h_l,h_l2,1:2:34,p_l,ge_l,e_l,end:end)),2),2.*std(squeeze(H(h_l,h_l2,1:2:34,p_l,ge_l,e_l,burn:end))'),...
%             marker{p_l},'MarkerSize',6,'MarkerFaceColor',colors{p_l})
        h(p_l).Color = colors{p_l}
        hold on 
    end 
    if e_l==1 && h_l==1 && h_l2==1
        title('HSD','FontWeight','normal','FontSize',FS)
    elseif e_l==2 && h_l==1 && h_l2==1
        title('HSG','FontWeight','normal','FontSize',FS)
    elseif e_l==3 && h_l==1 && h_l2==1
        title('CG','FontWeight','normal','FontSize',FS)
    end
if h_l==1 && h_l2==1
   ylim([0.4 1]) 
elseif h_l==1 && h_l2==2
   ylim([0 0.4]) 
elseif h_l==1 && h_l2==3
   ylim([0 0.4]) 
elseif h_l==2 && h_l2==1
   ylim([0 0.7]) 
elseif h_l==2 && h_l2==2
   ylim([0.0 1]) 
elseif h_l==2 && h_l2==3
   ylim([0 0.5]) 
end
% ylim([0 1])
if e_l==1
ylabel(strcat(h_str(h_l),'\rightarrow',h2_str(h_l2)))
end
yl=ylim
xlim([22 96])
xticks(30:10:90)
yticks(yl(1):0.2:yl(2))
set(gca,'FontName','Times New Roman','FontSize',FS);
end
set(gcf,'color','w')
if h_l==2 && h_l2==3
    I=legend('Protective','Detrimental','Harmful','Orientation','horizontal')
end
end
end
legend('boxoff')

%
I.FontSize=FS+1
newPosition = [0.45 0.001 0.1 0.07];
newUnits = 'normalized';
set(I,'Position', newPosition,'Units', newUnits);
grid off

% print('C:\Users\jbueren\Dropbox\habits\Draft\figures\transitions_all','-depsc')


%%
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

types_y=2
mean_wealth=reshape(mean_wealth{1},9,37,2,types_y,3);

%wealth distribution cond on good health
h_l=1
for p=1:3
colors = {  [0.4660    0.6740    0.1880]   [0.8500    0.3250    0.0980]   [0.8500    0.3250    0.0980] [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
pattern = { '-'  ':'  ':' '-.' '-'};
lw=[2.0 2.0]
FS=11
figure(2)
set(2,'position',[150    150    600    600])
for e_l=1:3
subplot(3,3,(p-1)*3+e_l)
for y_l=1:types_y
%     h(y_l)=plot(50:2:98,mean_wealth(3+p,13:37,y_l,e_l)./1000,'Color',colors{y_l},'linewidth',lw(y_l),'linestyle',pattern{y_l})
h(y_l)=plot(26:2:80,mean_wealth(4+p,1:28,h_l,y_l,e_l)./1000,'Color',colors{y_l},'linewidth',lw(y_l),'linestyle',pattern{y_l})
hold on
set(gca,'FontName','Times New Roman','FontSize',FS);
end
if p==1
    if e_l==1
    title('HSD',"FontWeight","normal")
%     ylim([0 150])
    elseif e_l==2
    title('HSG',"FontWeight","normal")
%     ylim([0 600])
    else
    title('CG',"FontWeight","normal")
%     ylim([0 1000])
    end
end

if p==1
    ylim([0 500])
elseif p==2
    ylim([0 1000])
else
    ylim([0 2000])
end

% ylim([0 2000])
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
print('C:\Users\jbueren\Dropbox\habits\Slides\v1\figures\wealth_moments','-depsc')
print('C:\Users\jbueren\Dropbox\habits\Draft\figures\wealth_moments','-depsc')


%% difference in median wealth across health states
for p=1:3
colors = {  [0.4660    0.6740    0.1880]   [0.9290    0.6940    0.1250]    [0.8500    0.3250    0.0980] [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
pattern = { '--'  '-'  ':' '-.' '-'};
lw=[1.5 1.5 2.5]
FS=11
figure(3)
set(3,'position',[950    150    750    750])
for e_l=1:3
subplot(3,3,(p-1)*3+e_l)
for y_l=1:3
h(y_l)=plot(26:2:80,(mean_wealth(4+p,1:28,1,y_l,e_l)-mean_wealth(4+p,1:28,2,y_l,e_l))./1000,'Color',colors{y_l},'linewidth',lw(y_l),'linestyle',pattern{y_l})
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
ylim([0 600])
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

end





%% Labor force participation by educ and health status

% clear all
close all
cd('C:\Users\jbueren\Google Drive\endo_health\metric_model\Results')
types_y=2
fileID=fopen('labor_force_participation.txt');
participation=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
participation=reshape(participation{1},5,37,2,3,types_y);
participation(5,:,1,1,1)
colors = {  [0.4660    0.6740    0.1880]   [0.9290    0.6940    0.1250]    [0.8500    0.3250    0.0980] [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
pattern = { '--'  '-'  ':' '-.' '-'};
lw=[1.5 1.5 2.5 1.5 2.5]
FS=11

figure(2)
set(2,'position',[150    150    750    350])
y_l=1

for e_l=1:3
subplot(1,3,e_l)
for h_l=1:2
h(h_l)=plot(25:2:65,participation(5,1:21,h_l,e_l,y_l),'Color',colors{h_l},'linewidth',lw(h_l),'linestyle',pattern{h_l})
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

%% Labor force participation by educ and health status & previous labor force

clear all
% close all
cd('C:\Users\jbueren\Google Drive\endo_health\metric_model\Results')
types_y=2
fileID=fopen('labor_force_participation_dynamic.txt');
participation=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
participation=reshape(participation{1},5,37,2,3,types_y,2);

colors = {  [0.4660    0.6740    0.1880]   [0.9290    0.6940    0.1250]    [0.8500    0.3250    0.0980] [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
pattern = { '--'  '-'  ':' '-.' '-'};
lw=[1.5 1.5 2.5 1.5 2.5]
FS=11

y_l=1
for f_l=1:2
figure(f_l+3)
set(f_l+3,'position',[150    150    750    350])
for e_l=1:3
subplot(1,3,e_l)
for h_l=1:2
h(h_l)=plot(26:2:60,participation(5,1:18,h_l,e_l,y_l,f_l),'Color',colors{h_l},'linewidth',lw(h_l),'linestyle',pattern{h_l})
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
I=legend('good health','bad health','Location','northwest','orientation','horizontal')
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
types_y=2
fileID=fopen('median_income.txt');
mean_wealth=textscan(fileID,'%14.10f','TreatAsEmpty',{'**************'});
fclose(fileID);
mean_wealth=reshape(mean_wealth{1},7,37,2,3,types_y);

colors = {  [0.4660    0.6740    0.1880]   [0.9290    0.6940    0.1250]    [0.8500    0.3250    0.0980] [0   0.4470    0.7410] [0.4940    0.1840    0.5560]};
pattern = { '--'  '-'  ':' '-.' '-'};
lw=[1.5 1.5 2.5]
FS=11
figure(5)
set(5,'position',[150    150    750    350])
y_l=1
for e_l=1:3
subplot(1,3,e_l)
for h_l=1:2
 plot(26:2:64,mean_wealth(7,1:20,h_l,e_l,y_l),'Color',colors{h_l},'linewidth',lw(h_l),'linestyle',pattern{h_l})
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
ylim([9.5 11.5])
xticks(0:10:60)
% yticks(0:25:150)
end
% I=legend('Protective','Detrimental','Harmful','Location','northwest','orientation','horizontal')
% legend('boxoff')
% I.FontSize=FS
% newPosition = [0.45 -0.02 0.1 0.1];
% newUnits = 'normalized';
% set(I,'Position', newPosition,'Units', newUnits);
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


