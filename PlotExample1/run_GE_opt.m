function [JOPT,JKF,JRDP_LMI,JRDP_PD,JIMM,JEM,nhyp_LMI,nhyp_PD,t_KF,t_LMI,t_PD,t_IMM,t_EM] = run_GE_opt(N,e)
%% Main code
opt = true;

%% Settings
tol = 1e-8;
%% Run Simulation
[x,y,sig,sys] = gil_ell_example1(N);
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
t_KF = zeros(N,1);

%RDP LMI
Pset_prev_LMI = 1:s;
nhyp_LMI = zeros(N+1,1);
% Time 0 RDP
for ii = 1:s
    xRDP_old_LMI(:,ii) = x0;
    P_RDP_old_LMI(:,:,ii) = P0;
    cRDP_old_LMI(ii) = ell0(ii);
end
nhyp_LMI(1) = s;
t_LMI = zeros(N,1);

%RDP PosDef (PD)
Pset_prev_PD = 1:s;
nhyp_PD = zeros(N+1,1);
% Time 0 RDP
for ii = 1:s
    xRDP_old_PD(:,ii) = x0;
    P_RDP_old_PD(:,:,ii) = P0;
    cRDP_old_PD(ii) = ell0(ii);
end
nhyp_PD(1) = s;
t_PD = zeros(N,1);


%IMM
mu_IMM = 1/s.*ones(s,1);
for ii = 1:s
    Phat_IMM(:,:,ii) = P_KF(:,:,1);
    xhat_IMM(:,ii) = xKF(:,1);
end
modeIMM = zeros(N+1,1);
modeIMM(1) = 1;
xIMM = zeros(n,N+1);
P_IMM = zeros(n,n,N+1);
% Time 0 IMM
xIMM(:,1) = xKF(:,1);
P_IMM(:,:,1) = P_KF(:,:,1);
t_IMM = zeros(N,1);

%% Run filters
% Time k 1,...,N
for k = 1:N
    % KF with perfect knowledge of mode
    tic
    [xKF(:,k+1),P_KF(:,:,k+1),cKF(k+1)] = MJLS_KF2(sys,xKF(:,k),P_KF(:,:,k),sig(k),sig(k+1),y(:,k+1),u(k),u(k+1),cKF(k),0);
    t_KF(k) = toc;
    
    % RDP using LMI
    tic
    [xRDP_LMI,P_RDP_LMI,cRDP_LMI,~,Pset_LMI] = MJLS_RDP_LMI_fast(sys,xRDP_old_LMI,P_RDP_old_LMI,cRDP_old_LMI,y(:,k+1),u(k),u(k+1),Pset_prev_LMI,e);
    t_LMI(k) = toc;
    xRDP_old_LMI = xRDP_LMI;
    P_RDP_old_LMI = P_RDP_LMI;
    cRDP_old_LMI = cRDP_LMI;
    Pset_prev_LMI = Pset_LMI;
    nhyp_LMI(k+1) = size(Pset_LMI,2);
    
    % RDP using PosDef
    tic
    [xRDP_PD,P_RDP_PD,cRDP_PD,~,Pset_PD] = MJLS_RDP_PosDefCheck(sys,xRDP_old_PD,P_RDP_old_PD,cRDP_old_PD,y(:,k+1),u(k),u(k+1),Pset_prev_PD,e);
    t_PD(k) = toc;
    xRDP_old_PD = xRDP_PD;
    P_RDP_old_PD = P_RDP_PD;
    cRDP_old_PD = cRDP_PD;
    Pset_prev_PD = Pset_PD;
    nhyp_PD(k+1) = size(Pset_PD,2);
    
    % IMM filter
    tic
    [xIMM(:,k+1),P_IMM(:,:,k+1),xhat_IMM,Phat_IMM,mu_IMM] = imm2(sys,xhat_IMM,Phat_IMM,mu_IMM,y(:,k+1),u(k),u(k+1));
    t_IMM(k) = toc;
    [~,modeIMM(k+1,1)] = max(mu_IMM);
end



%% KF using EM paper method
% EM
tic
[xEMs,sig_EM,~] = MJLS_EM(sys,y,u,sig,0,x);
t_EM = toc;

%% Optimal
if(opt)
    sig_pos = dec2bin(2^(N+1)-1:-1:0)-'0'+1;
    c = zeros(s^(N+1),1);
    for ii = 1:s^(N+1)
        sig_opt = sig_pos(ii,:);
        
        P_OPT = zeros(n,n,N+1);
        cOPT = zeros(N+1,1);
        xOPT = zeros(n,N+1);
        % Time 0 KF
        xOPT(:,1) = x0;
        P_OPT(:,:,1) = P0;
        cOPT(1) = ell0(sig_opt(1));
        
        for k = 1:N
            % KF
            [xOPT(:,k+1),P_OPT(:,:,k+1),cOPT(k+1)] = MJLS_KF2(sys,xOPT(:,k),P_OPT(:,:,k),sig_opt(k),sig_opt(k+1),y(:,k+1),u(k),u(k+1),cOPT(k));
        end
        c(ii,1) = cOPT(end);
    end
    [~,ind_opt] = min(c);
    sig_opt = sig_pos(ind_opt,:);
    
    P_OPT = zeros(n,n,N+1);
    cOPT = zeros(N+1,1);
    xOPT = zeros(n,N+1);
    % Time 0 KF
    xOPT(:,1) = x0;
    P_OPT(:,:,1) = P0;
    cOPT(1) = ell0(sig_opt(1));
    
    for k = 1:N
        % KF with optimal mode sequence
        [xOPT(:,k+1),P_OPT(:,:,k+1),cOPT(k+1)] = MJLS_KF2(sys,xOPT(:,k),P_OPT(:,:,k),sig_opt(k),sig_opt(k+1),y(:,k+1),u(k),u(k+1),cOPT(k));
    end
    xOPTs = MJLS_smooth2(sys,xOPT,P_OPT,u,sig_opt);
    JOPT = sum(cost_check2(sys,xOPTs,u,y,sig_opt));
    
    if(abs(JOPT-min(c))>tol || abs(JOPT - cOPT(end))>tol)
        disp('Error in OPT cost check');
        pause
    end
end

%% Smoothing
% RDP smooth procedure (calculate estimates using final Pset
[xRDP_post_PD,P_RDP_post_PD,cRDP_post_PD] = RDP_post(sys,Pset_PD(:,1),y,u);
[xRDP_post_LMI,P_RDP_post_LMI,cRDP_post_LMI] = RDP_post(sys,Pset_LMI(:,1),y,u);

% Smoothing
xKFs = MJLS_smooth2(sys,xKF,P_KF,u,sig);
xRDPs_LMI = MJLS_smooth2(sys,xRDP_post_LMI,P_RDP_post_LMI,u,Pset_LMI(:,1));
xRDPs_PD = MJLS_smooth2(sys,xRDP_post_PD,P_RDP_post_PD,u,Pset_PD(:,1));
xIMMs = MJLS_smooth2(sys,xIMM,P_IMM,u,modeIMM);

%% Cost checks
% Cost check
jkf = cost_check2(sys,xKFs,u,y,sig);
jrdp_lmi = cost_check2(sys,xRDPs_LMI,u,y,Pset_LMI(:,1));
jrdp_pd = cost_check2(sys,xRDPs_PD,u,y,Pset_PD(:,1));
jimm = cost_check2(sys,xIMMs,u,y,modeIMM);
jem = cost_check2(sys,xEMs,u,y,sig_EM);

% Costs
costKF = cKF(end);
costRDPs_LMI = cRDP_post_LMI(end);
costRDPs_PD = cRDP_post_PD(end);

JKF = sum(jkf);
JRDP_LMI = sum(jrdp_lmi);
JRDP_PD = sum(jrdp_pd);
JIMM = sum(jimm);
JEM = sum(jem);

nhyp_PD = sum(nhyp_PD);
nhyp_LMI = sum(nhyp_LMI);

t_KF = sum(t_KF);
t_LMI= sum(t_LMI);
t_PD = sum(t_PD);
t_IMM = sum(t_IMM);

if(abs(JKF-costKF)>tol)
    disp('Error in KF cost check');
    pause
end
if(abs(JRDP_LMI-costRDPs_LMI)>tol)
    disp('Error in RDP LMI cost check');
    pause
end
if(abs(JRDP_PD-costRDPs_PD)>tol)
    disp('Error in RDP PD cost check');
    pause
end


