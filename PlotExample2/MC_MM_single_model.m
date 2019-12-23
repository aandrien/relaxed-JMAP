%% Default code
clc; close all; clear;
addpath(genpath('../functions/'))
%% Main code
% Set noise generator for repeatability
% rng('shuffle');
rng(709587174,'twister') %seed used for paper
seed = rng;

%% Settings
MC = 500;
N  = 100;
e = 0.5;
T = 0.01;

%% Initialization
ne = length(e);

JKF = zeros(ne,1);
JRDP_PD = zeros(ne,1);
JEM = zeros(ne,1);
nHYP_PD = zeros(ne,1);

jkf = zeros(ne,MC);
jrdp_PD = zeros(ne,MC);
jem = zeros(ne,MC);
nhyp_PD = zeros(ne,MC);
t_KF = zeros(ne,MC);
t_PD = zeros(ne,MC);
t_EM = zeros(ne,MC);

%% Run MC simulations
for ii = 1:ne
    ii
    for jj = 1:MC
        fprintf('Iteration over e at %d \t of total %d \n \n Iteration over MC at %d \t of total %d \n \n',ii,ne,jj,MC)
        %         rng('shuffle');
%         [jkf(ii,jj),jrdp_PD(ii,jj),jem(ii,jj),nhyp_PD(ii,jj),t_KF(ii,jj),t_PD(ii,jj),t_EM(ii,jj)] = run_MM_single_model_errors(N,e(ii));
        tstart = tic;
        [rmsKF(:,jj),rmsRDP_PD(:,jj),rmsEM(:,jj),percRDP_PD(:,jj),percEM(:,jj),t_KF(:,jj),t_PD(:,jj),t_EM(:,jj)] = run_MM_single_model_errors(N,T,e(ii));
        tpassed = toc(tstart)
    end
end
save('MM_error_data3.mat')

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