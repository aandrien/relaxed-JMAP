function [x,y,u,sig,sys_mats] = example_MM_TT_full(N)
%Example of moving target with constant velocity and constant turn rates

%% Main code
%% Settings
% N = end time
dt = 0.1;
%% System
% Constant velocity
% A{1} = [1 dt 0 0;
%     0 1 0 0;
%     0 0 1 dt;
%     0 0 0 1];

% Constant turn rate
Omega = 0.1.*[0.1;20];
Omega(Omega==0) = []; %remove zero turn rate
nModes = length(Omega);
for ii = 1:nModes
    omega = Omega(ii);
    
    A{ii} = [1 sin(omega*dt)/omega 0 -(1-cos(omega*dt))/omega;
        0 cos(omega*dt) 0 -sin(omega*dt)
        0 (1-cos(omega*dt))/omega 1 sin(omega*dt)/omega
        0 sin(omega*dt) 0 cos(omega*dt)];
end

% Rest of system is same for each sigma
for jj = 1:nModes
    B{jj} = zeros(4,1);
    G{jj} = [0.5*dt^2 0; dt 0; 0 0.5*dt^2; 0 dt];
    %     G{jj} = eye(4);
    
    C{jj} = eye(4);
    D{jj} = 0;
    H{jj} = eye(4);
end

sys_mats.A = A; sys_mats.B = B; sys_mats.G = G;
sys_mats.C = C; sys_mats.D = D; sys_mats.H = H;

phi0 = 45;
v0 = 1;
x0 = [0;v0*cosd(phi0);0;v0*sind(phi0)];
% x0 = zeros(4,1);
P0 = diag([1^2,0.2^2,1^2,0.2^2]);
Q = 0.001.*eye(2);
R = 0.001.*eye(4);


rho0 = zeros(nModes,1);
for jj = 1:nModes
    rho0(jj) = 1/nModes;
end
% rho0 = [0.4; 0.3; 0.3];
rho0 = [0.9;0.1];
% rho_mat = [0.9 0.1 0 0 0;
%     0.1 0.8 0.1 0 0;
%     0 0.1 0.8 0.1 0;
%     0 0 0.1 0.8 0.1;
%     0 0 0 0.1 0.9];

% rho_mat = [0.4 0.3 0.3;
%     0.3 0.4 0.3;
%     0.3 0.3 0.4];

rho_mat = [0.6 0.3
    0.4 0.7];

% rho_mat = [0.85 0.1 0.3;
%     0.1 0.8 0.3;
%     0.05 0.1 0.4];


if(any(sum(rho_mat,1)~=1))
    disp('Error, probabilities do not sum up to one!')
    return
end

sys_mats.x0 = x0; sys_mats.P0 = P0; sys_mats.Q = Q; sys_mats.R = R;
sys_mats.rho0 = rho0; sys_mats.rho = rho_mat;

w = mvnrnd(zeros(1,2),Q,N);
% w = mvnrnd(zeros(1,4),Q,N);
v = mvnrnd(zeros(1,4),R,N);

%% Initialize variables
sig = zeros(N,1);
x = zeros(size(A{1},1),N);
y = zeros(size(C{1},1),N);
u = zeros(1,N);

%% Simulate model
for ii = 1:nModes
    sig(1+(ii-1)*round(N/nModes):ii*ceil(N/nModes)) = ii;
end

% sig(1) =  randsrc(1,1,[1,2,3;rho0(1),rho0(2),rho0(3)]); %Random variable with value 1 or 2 with prob of 0.5
x(:,1) = mvnrnd(x0,P0);
% sys_mats.x0 = x(1,:);
for k = 1:N
    if(k<N)
%         x(:,k+1) = A{sig(k)}*x(:,k)+B{sig(k)}*u(:,k)+G{sig(k)}*w(k,:)';
                x(:,k+1) = A{sig(k)}*x(:,k)+B{sig(k)}*u(:,k);
        
%         if(sig(k) == 1)
%             sig(k+1) = randsrc(1,1,[1,2,3;rho_mat(1,1),rho_mat(2,1),rho_mat(3,1)]); %Random variable 1 with prob 0.6 and 2 with prob 0.4
%         elseif(sig(k) == 2)
%             sig(k+1) = randsrc(1,1,[1,2,3;rho_mat(1,2),rho_mat(2,2),rho_mat(3,2)]); %Random variable 1 with prob 0.15 and 2 with prob 0.85
%         elseif(sig(k) == 3)
%             sig(k+1) = randsrc(1,1,[1,2,3;rho_mat(1,3),rho_mat(2,3),rho_mat(3,3)]); %Random variable 1 with prob 0.15 and 2 with prob 0.85
%         end
    end
    y(:,k) = C{sig(k)}*x(:,k)+H{sig(k)}*v(k,:)';
end
