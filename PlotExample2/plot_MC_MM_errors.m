%% Default code
clc; close all; clear;



%% Main code
load('MM_error_data2.mat')
varKF = var(rmsKF,[],2)
varEM = var(rmsEM,[],2)
varRDP = var(rmsRDP_PD,[],2)

rmsKF = mean(rmsKF,2)
rmsEM = mean(rmsEM,2)
rmsRDP = mean(rmsRDP_PD,2)

percRDP = mean(percRDP_PD)
percEM = mean(percEM)

%% Figures
clearvars;
saveplot = true;
%% Figure Settings
lineStyles = linspecer(4);

load('MM_plot_instance2.mat')

figure
plot(x(1,:),x(3,:),'color',lineStyles(1,:))
hold on
plot(xKFs(1,:),xKFs(3,:),'color',lineStyles(2,:))
plot(xRDPs(1,:),xRDPs(3,:),'color',lineStyles(3,:))
% plot(e,JIMM./JOPT)
plot(xEMs(1,:),xEMs(3,:),'color',lineStyles(4,:))
ylabel('ylabel')
xlabel('xlabel')
legend('State','KF','PosDef','EM','location','SouthEast')
grid on

if(saveplot)
    set(gcf,'color','w');
    print -depsc ..\..\..\V9\Pictures\MM_instance.eps
end
