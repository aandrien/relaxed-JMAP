function [x,y,sig,sys_mats] = gil_ell_example2D(N)
%Example 1 from 'An online sequential algorithm for the estimation of
%transition probabilities for jump Markov linear systems 2006 Orguner'

%% Main code
%% Settings
% N = end time

%% System
A{1} = eye(2); A{2} = eye(2); 
B{1} = zeros(2,1); B{2} = zeros(2,1);
G{1} = eye(2); G{2} = eye(2);
C{1} = zeros(2,2); C{2} = eye(2);
D{1} = zeros(2,1); D{2} = zeros(2,1);
H{1} = 100.*eye(2); H{2} = 10.*eye(2);
% H{1} = 10; H{2} = 10;
sys_mats.A = A; sys_mats.B = B; sys_mats.G = G;
sys_mats.C = C; sys_mats.D = D; sys_mats.H = H;

x0 = [100;80];
P0 = diag([20^2,20^2]);
Q = diag([2^2,2^2]);
R = diag([1,1]);

rho0(1) = 0.2;
rho0(2) = 0.8;

% rho_mat = [0.6, 0.15
%     0.4, 0.85];

rho_mat = [0.15, 0.3
    0.85, 0.7];


sys_mats.x0 = x0; sys_mats.P0 = P0; sys_mats.Q = Q; sys_mats.R = R;
sys_mats.rho0 = rho0; sys_mats.rho = rho_mat;

w = mvnrnd(zeros(1,2),Q,N);
v = mvnrnd(zeros(1,2),R,N);

%% Initialize variables
sig = zeros(N,1);
x = zeros(size(A{1},1),N);
y = zeros(size(C{1},1),N);

%% Simulate model
sig(1) =  randsrc(1,1,[1,2;rho0(1),rho0(2)]); %Random variable with value 1 or 2 with prob of 0.5
x(:,1) = mvnrnd(x0,P0);
% sys_mats.x0 = x(1);

for k = 1:N
    if(k<N)
        x(:,k+1) = A{sig(k)}*x(:,k)+G{sig(k)}*w(k,:)';
        
        if(sig(k) == 1)
            sig(k+1) = randsrc(1,1,[1,2;rho_mat(1,1),rho_mat(2,1)]); %Random variable 1 with prob 0.6 and 2 with prob 0.4
%             sig(k+1) = 1;
        elseif(sig(k) == 2)
            sig(k+1) = randsrc(1,1,[1,2;rho_mat(1,2),rho_mat(2,2)]); %Random variable 1 with prob 0.15 and 2 with prob 0.85
%             sig(k+1) = 1;
        end
    end
    y(:,k) = C{sig(k)}*x(:,k)+H{sig(k)}*v(k,:)';
end