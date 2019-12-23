%% Default code
clc; close all; clear;

saveplot = true;

%% Figure Settings
% figure_configuration_IEEE_standard
set(0, 'DefaultLineLineWidth', 2);
set(0,'DefaultAxesXGrid','on')
set(0,'DefaultAxesYGrid','on')
set(0,'defaultAxesFontName','Times New Roman');
set(0,'defaultAxesFontSize',12);
set(0,'defaultTextFontName','Times New Roman');
set(0,'defaultTextFontSize',12);
set(0,'defaultLegendFontName','Times New Roman');
set(0,'defaultLegendFontSize',12);

lineStyles = linspecer(4);

%% Load Data
load('GE_opt_data3.mat')

JOPT = mean(jopt,2);
JKF = mean(jkf,2);
JRDP_LMI = mean(jrdp_LMI,2);
JRDP_PD = mean(jrdp_PD,2);
JIMM = mean(jimm,2);
JEM = mean(jem,2);
nHYP_PD = mean(nhyp_PD,2);
nHYP_LMI = mean(nhyp_LMI,2);

T_KF = mean(t_KF,2);
T_LMI = mean(t_LMI,2);
T_PD = mean(t_PD,2);
T_EM = mean(t_EM,2);
T_IMM = mean(t_IMM,2);

%% Figures
%Figure 1
figure
ax(1) = subplot(2,1,1);
plot(e,JKF./JOPT,'color',lineStyles(1,:))
hold on
plot(e,JRDP_LMI./JOPT,'color',lineStyles(2,:))
plot(e,JRDP_PD./JOPT,'color',lineStyles(3,:),'LineStyle','--')
% plot(e,JIMM./JOPT)
plot(e,JEM./JOPT,'color',lineStyles(4,:))
ylim([0.9,1.6])
ylabel('Scaled cost [-]','interpreter','latex')
legend('KF','LMI','PosDef','EM')
grid on
ax(2) = subplot(2,1,2);
plot(e,nHYP_LMI,'color',lineStyles(2,:))
hold on
plot(e,nHYP_PD,'color',lineStyles(3,:),'LineStyle','--')
grid on
legend('LMI','PosDef')
ylabel('Number of hypotheses [-]','interpreter','latex')



xlabel('Relaxation parameter $\epsilon$','interpreter','latex')
linkaxes(ax,'x')
xlim([0, 5])

if(saveplot)
    set(gcf,'color','w');
    export_fig ..\..\..\V8\Pictures\MC_cost_vs_epsilon.eps
end

figure
semilogy(e,T_KF,'color',lineStyles(1,:))
hold on
semilogy(e,T_LMI,'color',lineStyles(2,:))
semilogy(e,T_PD,'color',lineStyles(3,:))
semilogy(e,T_EM,'color',lineStyles(4,:))
grid on
legend('KF','LMI','PosDef','EM')
ylabel('Computation time [s]','interpreter','latex')
xlabel('Relaxation parameter $\epsilon$','interpreter','latex')
xlim([0, 1])

if(saveplot)
    set(gcf,'color','w');
    export_fig ..\..\..\V8\Pictures\MC_time_vs_epsilon.eps
end





%% Figures

% figure(1)
% subplot(3,1,1)
% plot(t(plot_t),p_MPC((plot_t),2),'color',lineStyles(1,:))
% hold on
% plot(t(plot_t),p_base((plot_t),2),'color',lineStyles(2,:))
% plot(t(plot_t),p_ref((plot_t),3),'color',lineStyles(3,:))
% 
% legend('MPC','Saturated PD + FF','Reference')
% ylabel('Position $[m]$','interpreter','latex')
% 
% subplot(3,1,2)
% plot(t(plot_t),p_MPC((plot_t),2)-p_ref((plot_t),3),'color',lineStyles(1,:))
% hold on
% plot(t(plot_t),p_base((plot_t),2)-p_ref((plot_t),3),'color',lineStyles(2,:))
% ylabel('Position Error $[m]$','interpreter','latex')
% 
% subplot(3,1,3)
% plot(t(plot_t),a_MPC((plot_t),2),'color',lineStyles(1,:))
% hold on
% plot(t(plot_t),a_base((plot_t),2),'color',lineStyles(2,:))
% plot([t(1) t(plot_t(end))],[L L],'k--')
% plot([t(1) t(plot_t(end))],[-L -L],'k--')
% ylim([-1.1*L 1.1*L])
% ylabel('Acceleration $[\frac{m}{s^2}]$','interpreter','latex')
% xlabel('Time $[s]$','interpreter','latex')
% 
% if(saveplot)
%     set(gcf,'color','w');
%     export_fig ..\..\V8\Pictures\MPCvPD_sine.eps
% end
