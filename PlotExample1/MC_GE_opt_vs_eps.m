%% Default code
clc; close all; clear;
addpath(genpath('../functions/'))
%% Main code
% Set noise generator for repeatability
rng('shuffle');

rng(701674144,'twister') %seed used for plot in paper
seed = rng;

%% Settings
MC = 50;
N  = 12;
e = [0:0.1:1,2:1:10];
% algos = {'RDP','KF','RDP_LMI','OPT','IMM','EM'};
algos = {'RDP','RDP_LMI','OPT','KF'};

%% Initialization
ne = length(e);

JOPT = zeros(ne,1);
JKF = zeros(ne,1);
JRDP_PD = zeros(ne,1);
JRDP_LMI = zeros(ne,1);
JIMM = zeros(ne,1);
JEM = zeros(ne,1);
nHYP_PD = zeros(ne,1);
nHYP_LMI = zeros(ne,1);

jopt = zeros(ne,MC);
jkf = zeros(ne,MC);
jrdp_PD = zeros(ne,MC);
jrdp_LMI = zeros(ne,MC);
jimm = zeros(ne,MC);
jem = zeros(ne,MC);
nhyp_PD = zeros(ne,MC);
nhyp_LMI = zeros(ne,MC);
t_KF = zeros(ne,MC);
t_LMI = zeros(ne,MC);
t_PD = zeros(ne,MC);
t_IMM = zeros(ne,MC);
t_EM = zeros(ne,MC);

%% Run MC simulations
for ii = 1:ne
    ii
    for jj = 1:MC
        fprintf('Iteration over e at %d \t of total %d \n \n Iteration over MC at %d \t of total %d \n \n',ii,ne,jj,MC)
        %         rng('shuffle');
        [jopt(ii,jj),jkf(ii,jj),jrdp_LMI(ii,jj),jrdp_PD(ii,jj),jimm(ii,jj),jem(ii,jj),nhyp_LMI(ii,jj),nhyp_PD(ii,jj),t_KF(ii,jj),t_LMI(ii,jj),t_PD(ii,jj),t_IMM(ii,jj),t_EM(ii,jj)] = run_GE_opt(N,e(ii));
    end
%     JOPT(ii) = sum(jopt)/MC;
%     JKF(ii) = sum(jkf)/MC;
%     JRDP_PD(ii) = sum(jrdp_PD)/MC;
%     JRDP_LMI(ii) = sum(jrdp_LMI)/MC;
%     JIMM(ii) = sum(jimm)/MC;
%     JEM(ii) = sum(jem)/MC;
%     nHYP_PD(ii) = sum(nhyp_PD)/MC;
%     nHYP_LMI(ii) = sum(nhyp_LMI)/MC;
end

if(~any(strcmp(algos,'OPT')))
    JOPT= ones(ne,1);
end

save('GE_opt_data3.mat')
%% Figures
figure
subplot(2,1,1)
plot(e,JKF./JOPT)
hold on
plot(e,JRDP_PD./JOPT)
plot(e,JRDP_LMI./JOPT)
plot(e,JIMM./JOPT)
plot(e,JEM./JOPT)
ylabel('$\frac{J}{J^*}$','interpreter','latex','Fontsize',20,'Rotation',0)
legend('KF','RDP','RDP_{LMI}','IMM','EM')
grid on
subplot(2,1,2)
plot(e,nHYP_PD)
hold on
plot(e,nHYP_LMI)
grid on
legend('RDP','RDP_{LMI}')
ylabel('Average number of hypotheses over horizon')
xlabel('$\epsilon$','interpreter','latex','Fontsize',20)