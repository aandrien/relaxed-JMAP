function [x,y,u,sig,sys_mats] = example_MM_TT(N)
%Example of moving target with constant velocity and constant turn rates

%% Main code
%% Settings
% N = end time
dt = 0.1;
%% System
% Constant velocity
A{1} = [1 dt 0 0;
    0 1 0 0;
    0 0 1 dt;
    0 0 0 1];

% Constant turn rate
Omega = 1.*[-2;2];
Omega(Omega==0) = []; %remove zero turn rate
nModes = length(Omega)+1;
for ii = 2:nModes
    omega = Omega(ii-1);
    
    A{ii} = [1 sin(omega*dt)/omega 0 -(1-cos(omega*dt))/omega;
        0 cos(omega*dt) 0 -sin(omega*dt)
        0 (1-cos(omega*dt))/omega 1 sin(omega*dt)/omega
        0 sin(omega*dt) 0 cos(omega*dt)];
end

% Rest of system is same for each sigma
for jj = 1:nModes
    B{jj} = zeros(4,1);
    %     G{jj} = [0.5*dt^2 0; dt 0; 0 0.5*dt^2; 0 dt];
    G{jj} = eye(4);
    
    C{jj} = [1 0 0 0;
        0 0 1 0];
    D{jj} = 0;
    H{jj} = eye(2);
end

sys_mats.A = A; sys_mats.B = B; sys_mats.G = G;
sys_mats.C = C; sys_mats.D = D; sys_mats.H = H;

phi0 = 45;
v0 = 1;
x0 = [0;v0*cosd(phi0);0;v0*sind(phi0)];
% x0 = zeros(4,1);
P0 = diag([1^2,0.2^2,1^2,0.2^2]);
Q = 0.001.*eye(4);
R = 0.01.*eye(2);


rho0 = zeros(nModes,1);
for jj = 1:nModes
    rho0(jj) = 1/nModes;
end
rho0 = [0.6; 0.2; 0.2];
% rho_mat = [0.9 0.1 0 0 0;
%     0.1 0.8 0.1 0 0;
%     0 0.1 0.8 0.1 0;
%     0 0 0.1 0.8 0.1;
%     0 0 0 0.1 0.9];

% rho_mat = [0.4 0.3 0.3;
%     0.3 0.4 0.3;
%     0.3 0.3 0.4];

% rho_mat = [0.7 0.15 0.05;
%     0.25 0.7 0.25;
%     0.05 0.15 0.7];

rho = [0.85 0.1 0.05;
    0.1 0.8 0.1;
    0.05 0.1 0.85];


if(any(sum(rho,1)~=1))
    disp('Error, probabilities do not sum up to one!')
    return
end

ell0 = -log(rho0);
ell = -log(rho);
ell(ell==Inf) = 5e1; %replace INF with high number, will cause numerical problems otherwise

%% Add variables to system structure
sys_mats.x0 = x0; sys_mats.P0 = P0; sys_mats.Q = Q; sys_mats.R = R;
sys_mats.ell0 = ell0; sys_mats.ell = ell;
sys_mats.rho0 = rho0; sys_mats.rho = rho;

sys_mats.n = size(A{1},1);
sys_mats.s = length(A);
sys_mats.N = N;
sys_mats.r = size(C{1},1);

% w = mvnrnd(zeros(1,2),Q,N);
w = mvnrnd(zeros(1,4),Q,N);
v = mvnrnd(zeros(1,2),R,N);

%% Initialize variables
sig = zeros(N+1,1);
x = zeros(size(A{1},1),N);
y = zeros(size(C{1},1),N);
u = zeros(1,N+1);

%% Simulate model
% for ii = 1:nModes
%     sig(1+(ii-1)*round(N/nModes):ii*ceil(N/nModes)) = ii;
% end

sig(1) =  randsrc(1,1,[1,2,3;rho0(1),rho0(2),rho0(3)]); %Random variable with value 1 or 2 with prob of 0.5
x(:,1) = mvnrnd(x0,P0);
% sys_mats.x0 = x(1,:);
for k = 2:N+1
    if(sig(k-1) == 1)
        sig(k) = randsrc(1,1,[1,2,3;rho(1,1),rho(2,1),rho(3,1)]); %Random variable 1 with prob 0.6 and 2 with prob 0.4
    elseif(sig(k-1) == 2)
        sig(k) = randsrc(1,1,[1,2,3;rho(1,2),rho(2,2),rho(3,2)]); %Random variable 1 with prob 0.15 and 2 with prob 0.85
    elseif(sig(k-1) == 3)
        sig(k) = randsrc(1,1,[1,2,3;rho(1,3),rho(2,3),rho(3,3)]); %Random variable 1 with prob 0.15 and 2 with prob 0.85
    end
        x(:,k) = A{sig(k)}*x(:,k-1)+B{sig(k)}*u(:,k-1)+G{sig(k)}*w(k-1,:)';
%     x(:,k) = A{sig(k)}*x(:,k-1)+B{sig(k)}*u(:,k-1);

    y(:,k) = C{sig(k)}*x(:,k)+H{sig(k)}*v(k-1,:)';
end
