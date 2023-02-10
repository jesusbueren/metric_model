clear all
close all
cd("C:\Users\jbueren\Google Drive\endo_health\metric_model\Results")
FS=9
fileID = fopen('montecarlo_beta.txt','r');
beta = fscanf(fileID,'%f')
fclose(fileID);
fileID = fopen('montecarlo_nu.txt','r');
nu = fscanf(fileID,'%f')
fclose(fileID);
fileID = fopen('montecarlo_w.txt','r');
w = fscanf(fileID,'%f')
fclose(fileID);
fileID = fopen('montecarlo_0.txt','r');
ini = fscanf(fileID,'%f')
fclose(fileID);
fileID = fopen('montecarlo_rho.txt','r');
rho = fscanf(fileID,'%f')
fclose(fileID);

beta_true=8.84794
nu_true=0.02
w_true=0.01
ini_true=0.25  
rho_true=0.9

for f_l=1:5
    if f_l==1
        h=beta
        t=beta_true
    elseif f_l==2
        h=nu
        t=nu_true
    elseif f_l==3
        h=w
        t=w_true
    elseif f_l==4
        h=ini
        t=ini_true
    elseif f_l==5
    h=rho
    t=rho_true
    end
        
    subplot(2,3,f_l)
    histogram(h)
    yl = ylim()
    hold on
    line([t, t], [yl(1), yl(2)],'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 2);
     if f_l==1
        title('\beta')
    elseif f_l==2
        title('\nu')
    elseif f_l==3
        title('w')
    elseif f_l==4
        title('0')
    elseif f_l==5
        title('\rho')
    end
    
end
