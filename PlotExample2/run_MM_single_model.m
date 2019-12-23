function [JKF,JRDP_PD,JEM,nhyp_PD,t_KF,t_PD,t_EM] = run_MM_single_model(N,e)
%% Main code
opt = true;

%% Settings
tol = 1e-8;
%% Run Simulation
[x,y,u,sig,sys] = example_MM_TT(N);
% u = zeros(1,N+1);

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


%% Run filters
% Time k 1,...,N
for k = 1:N
    % KF with perfect knowledge of mode
    tic
    [xKF(:,k+1),P_KF(:,:,k+1),cKF(k+1)] = MJLS_KF2(sys,xKF(:,k),P_KF(:,:,k),sig(k),sig(k+1),y(:,k+1),u(k),u(k+1),cKF(k),0);
    t_KF(k) = toc;
    
    % RDP using PosDef
    tic
    [xRDP_PD,P_RDP_PD,cRDP_PD,~,Pset_PD] = MJLS_RDP_PosDefCheck(sys,xRDP_old_PD,P_RDP_old_PD,cRDP_old_PD,y(:,k+1),u(k),u(k+1),Pset_prev_PD,e);
    t_PD(k) = toc;
    xRDP_old_PD = xRDP_PD;
    P_RDP_old_PD = P_RDP_PD;
    cRDP_old_PD = cRDP_PD;
    Pset_prev_PD = Pset_PD;
    nhyp_PD(k+1) = size(Pset_PD,2);
    
end



%% KF using EM paper method
% EM
tic
[xEMs,sig_EM,~] = MJLS_EM(sys,y,u,sig,0,x);
t_EM = toc;

%% Smoothing
% RDP smooth procedure (calculate estimates using final Pset
[xRDP_post_PD,P_RDP_post_PD,cRDP_post_PD] = RDP_post(sys,Pset_PD(:,1),y,u);

% Smoothing
xKFs = MJLS_smooth2(sys,xKF,P_KF,u,sig);
xRDPs_PD = MJLS_smooth2(sys,xRDP_post_PD,P_RDP_post_PD,u,Pset_PD(:,1));

%% Cost checks
% Cost check
jkf = cost_check2(sys,xKFs,u,y,sig);
jrdp_pd = cost_check2(sys,xRDPs_PD,u,y,Pset_PD(:,1));
jem = cost_check2(sys,xEMs,u,y,sig_EM);

% Costs
costKF = cKF(end);
costRDPs_PD = cRDP_post_PD(end);

JKF = sum(jkf);
JRDP_PD = sum(jrdp_pd);
JEM = sum(jem);

nhyp_PD = sum(nhyp_PD);

t_KF = sum(t_KF);
t_PD = sum(t_PD);

if(abs(JKF-costKF)>tol)
    disp('Error in KF cost check');
    pause
end
if(abs(JRDP_PD-costRDPs_PD)>tol)
    disp('Error in RDP PD cost check');
    pause
end


