%% Default code
clc; close all; clear all;
addpath(genpath('../functions/'))
%% Main code
opt = false;
%% Settings
N = 100; % end time
e =0.5;
T = 0.01;

% Set noise generator for repeatability
% rng('shuffle');
rng(709587174,'twister') %seed used for paper
seed = rng;


%% Run Simulation
% [x,y,u,sig,sys] = example_MM_TT(N);
[x,y,u,sig,sys] = example_MM_TT_scaled(N,T);
u = zeros(1,N+1);

n = sys.n;
s = sys.s;
x0 = sys.x0;
P0 = sys.P0;
ell0 = sys.ell0;

%% Initialize filter
%KF
P_KF = zeros(n,n,N+1);
P_KF(:,:,1) = eye(n);
cKF = zeros(N+1,1);
xKF = zeros(n,N+1);
% Time 0 KF
xKF(:,1) = x0;
P_KF(:,:,1) = P0;
cKF(1) = ell0(sig(1));

%RDP
Pset_prev = 1:s;
xRDP_old = x0;
P_RDP_old(:,:,1) = P0;
cRDP_old = zeros(s,1);
modeRDP = zeros(N+1,1);
xhat_RDP = zeros(n,N+1);
nhyp = zeros(N+1,1);
% Time 0 RDP
for ii = 1:s
    xRDP_old(:,ii) = x0;
    P_RDP_old(:,:,ii) = P0;
    cRDP_old(ii) = ell0(ii);
end
[~,ind] = min(cRDP_old);
modeRDP(1) = ind;
xhat_RDP(:,1) = xRDP_old(:,ind);
Phat_RDP(:,:,1) = P_RDP_old(:,:,ind);
cRDP_hat(1,1) = cRDP_old(ind);

%IMM
mu_IMM = 1/s.*ones(s,1);
for ii = 1:s
    Phat_IMM(:,:,ii) = P_KF(:,:,1);
    xhat_IMM(:,ii) = xKF(:,1);
end
modeIMM(1) = 1;
xIMM = zeros(n,N+1);
P_IMM = zeros(n,n,N+1);
% Time 0 IMM
xIMM(:,1) = xKF(:,1);
P_IMM(:,:,1) = P_KF(:,:,1);

%% Run filters
% Time k 1,...,N
for k = 1:N
    % KF with perfect knowledge of mode
    [xKF(:,k+1),P_KF(:,:,k+1),cKF(k+1)] = MJLS_KF2(sys,xKF(:,k),P_KF(:,:,k),sig(k),sig(k+1),y(:,k+1),u(k),u(k+1),cKF(k));
    
    % RDP filter procedure (keeping the best at time k)
    tic
%     [xRDP,P_RDP,cRDP,modeRDP(k+1),Pset] = MJLS_RDP_minmax2(xRDP_old,P_RDP_old,cRDP_old,y(:,k+1),u(k),u(k+1),Pset_prev,e,eta);
%     [xRDP,P_RDP,cRDP,modeRDP(k+1),Pset] = MJLS_RDP_LMI_fast(sys,xRDP_old,P_RDP_old,cRDP_old,y(:,k+1),u(k),u(k+1),Pset_prev,e);
%     [xRDP,P_RDP,cRDP,modeRDP(k+1),Pset] = MJLS_RDP_LMI(xRDP_old,P_RDP_old,cRDP_old,y(:,k+1),u(k),u(k+1),Pset_prev,e);
    [xRDP,P_RDP,cRDP,modeRDP(k+1),Pset] = MJLS_RDP_PosDefCheck(sys,xRDP_old,P_RDP_old,cRDP_old,y(:,k+1),u(k),u(k+1),Pset_prev,e);
    t(k) = toc;
    xRDP_old = xRDP;
    P_RDP_old = P_RDP;
    cRDP_old = cRDP;
    Pset_prev = Pset;
    xhat_RDP(:,k+1) = xRDP(:,1);
    Phat_RDP(:,:,k+1) = P_RDP(:,:,1);
    cRDP_hat(k+1,1) = cRDP(1);
    size(Pset)
    nhyp(k+1) = size(Pset,2);
    
    % IMM filter
    [xIMM(:,k+1),P_IMM(:,:,k+1),xhat_IMM,Phat_IMM,mu_IMM] = imm2(sys,xhat_IMM,Phat_IMM,mu_IMM,y(:,k+1),u(k),u(k+1));
    [~,modeIMM(k+1,1)] = max(mu_IMM);
end



%% KF using EM paper method
% EM
[xEMs,sig_EM,nEM] = MJLS_EM(sys,y,u,Pset(:,1),0,x);

%% Optimal
if(opt)
    %     sig_pos = dec2bin(2^N-1:-1:0)-'0'+1;
    sig_pos = dec2base(0:(3^N-1),3)-48+1;
    tic
    for ii = 1:s^N
        sig_opt = sig_pos(ii,:);
        
        P_OPT = zeros(n,n,N);
        cOPT = zeros(N,1);
        xOPT = zeros(n,N);
        % Time 0 KF
        C0 = C{sig_opt(1)};
        D0 = D{sig_opt(1)};
        H0 = H{sig_opt(1)};
        R0 = H0*R*H0';
        L0 = inv(R0+C0*P0*C0');
        K0 = P0*C0'*L0;
        xOPT(:,1) = x0+K0*(y(:,1)-C0*x0-D0*u(:,1));
        P_OPT(:,:,1) = P0-K0*C0*P0;
        cOPT(1) = 0.5*mnorm(y(:,1)-C0*x0-D0*u(:,1),L0)^2+ell0(sig_opt(1));
        
        for k = 1:N-1
            % KF
            [xOPT(:,k+1),P_OPT(:,:,k+1),cOPT(k+1)] = MJLS_KF(xOPT(:,k),P_OPT(:,:,k),sig_opt(k),sig_opt(k+1),y(:,k+1),u(k),u(k+1),cOPT(k));
        end
        c(ii) = cOPT(end);
    end
    avg_time_opt = toc/ii
    [~,ind_opt] = min(c);
    sig_opt = sig_pos(ind_opt,:);
    
    P_OPT = zeros(n,n,N);
    cOPT = zeros(N,1);
    xOPT = zeros(n,N);
    % Time 0 KF
    C0 = C{sig_opt(1)};
    D0 = D{sig_opt(1)};
    H0 = H{sig_opt(1)};
    R0 = H0*R*H0';
    L0 = inv(R0+C0*P0*C0');
    K0 = P0*C0'*L0;
    xOPT(:,1) = x0+K0*(y(:,1)-C0*x0-D0*u(:,1));
    P_OPT(:,:,1) = P0-K0*C0*P0;
    cOPT(1) = 0.5*mnorm(y(:,1)-C0*x0-D0*u(:,1),L0)^2+ell0(sig_opt(1));
    
    for k = 1:N-1
        % KF with perfect knowledge of mode
        [xOPT(:,k+1),P_OPT(:,:,k+1),cOPT(k+1)] = MJLS_KF(xOPT(:,k),P_OPT(:,:,k),sig_opt(k),sig_opt(k+1),y(:,k+1),u(k),u(k+1),cOPT(k));
    end
    xOPTs = MJLS_smooth(xOPT,P_OPT,u,sig_opt);
    JOPT = cost_check(xOPTs,u,y,sig_opt);
    eOPT = xOPT-x;
    eOPTs = xOPTs-x;
    eOPTrms = rms(eOPT)
    eOPTrms_s = rms(eOPTs)
    mae_OPT = mae_vec(xOPT,x)
    mae_OPTs = mae_vec(xOPTs,x)
    costOPT= min(c)
end


%% Smoothing
% RDP smooth procedure (calculate estimates using final Pset
[xRDP_post,P_RDP_post,cRDP_post] = RDP_post(sys,Pset(:,1),y,u);
[xhat_RDP,Phat_RDP,cRDP_hat] = RDP_post(sys,modeRDP,y,u);

% Smoothing
xKFs = MJLS_smooth2(sys,xKF,P_KF,u,sig);
xRDPs = MJLS_smooth2(sys,xRDP_post,P_RDP_post,u,Pset(:,1));
xIMMs = MJLS_smooth2(sys,xIMM,P_IMM,u,modeIMM);

plot(x(1,:),x(3,:))
hold on
% plot(y(1,:),y(2,:))
plot(xKFs(1,:),xKFs(3,:))
plot(xEMs(1,:),xEMs(3,:))
plot(xIMMs(1,:),xIMMs(3,:))
plot(xRDPs(1,:),xRDPs(3,:))
legend('x','KFs','EM','IMM','RDP')

%% Cost checks
% Cost check
JKF = cost_check2(sys,xKFs,u,y,sig);
JRDP = cost_check2(sys,xRDPs,u,y,Pset);
JIMM = cost_check2(sys,xIMMs,u,y,modeIMM);
JEM = cost_check2(sys,xEMs,u,y,sig_EM);

%% Rescale
xKFs = T.*xKFs;
xRDPs = T.*xRDPs;
xEMs = T.*xEMs;
x = T.*x;

% Errors
eRDP= xhat_RDP-x;
eRDPs= xRDPs-x;
eIMM= xIMM-x;
eIMMs= xIMMs-x;
eKF= xKF-x;
eKFs = xKFs-x;

eRDPrms = rms(eRDP);
eRDPrms_s = rms(eRDPs);
eIMMrms = rms(eIMM);
eIMMrms_s = rms(eIMMs);
eKFrms = rms(eKF);
eKFrms_s = rms(eKFs);


mae_RDP = mae_vec(xhat_RDP,x);
mae_RDPs = mae_vec(xRDPs,x);
mae_IMM = mae_vec(xIMM,x);
mae_IMMs = mae_vec(xIMMs,x);
mae_KF = mae_vec(xKF,x);
mae_KFs = mae_vec(xKFs,x);

costKF = cKF(end)
costRDP = cRDP_hat(end)
costRDPs = cRDP_post(end)
costIMMs = sum(JIMM)
costEMs = sum(JEM)

% sig = sig(1:N);
percRDP_filt = sum(modeRDP ==sig)/size(sig,1)
percRDP_smooth = sum(Pset(:,1) ==sig)/size(sig,1)
percIMM = sum(modeIMM==sig)/size(sig,1)
percEMM = sum(sig_EM==sig)/size(sig,1)

% save('MM_plot_instance3.mat','xKFs','xEMs','xRDPs','x')
