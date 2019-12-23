function [x,y,u,sig,sys_mats] = example_MM_track(N)
%Example 2 from 'An online sequential algorithm for the estimation of
%transition probabilities for jump Markov linear systems 2006 Orguner'

%% Main code
%% Settings
% N = end time
dt = 0.1;
%% System
A{1} = [1 dt 0 0;
    0 1 0 0;
    0 0 1 dt;
    0 0 0 1]; 
A{2} = A{1}; A{3} = A{1};
B{1} = zeros(4,1);
B{2} = [-1.225;-0.35;1.225;0.35]; B{3} = -B{2};
G{1} = eye(4); G{2} = eye(4); G{3} = eye(4);
C{1} = eye(4); C{2} = eye(4); C{3} = eye(4);
D{1} = 0; D{2} = 0; D{3} = 0;
H{1} = eye(4); H{2} = eye(4); H{3} = eye(4);
% H{1} = 10; H{2} = 10;
sys_mats.A = A; sys_mats.B = B; sys_mats.G = G;
sys_mats.C = C; sys_mats.D = D; sys_mats.H = H;

x0 = zeros(4,1);
P0 = 2^2*eye(4);
Q = 00001.*0.1^2.*eye(4);
R = 100.*eye(4);

rho0(1) = 0.5;
rho0(2) = 0.25;
rho0(3) = 0.25;

rho_mat = [0.9 0.05 0.05;
    0.05 0.9 0.05;
    0.05 0.05 0.9];

if(any(sum(rho_mat,1)~=1))
    disp('Error, probabilities do not sum up to one!')
    return
end

sys_mats.x0 = x0; sys_mats.P0 = P0; sys_mats.Q = Q; sys_mats.R = R;
sys_mats.rho0 = rho0; sys_mats.rho = rho_mat;

w = mvnrnd(zeros(1,4),Q,N);
v = mvnrnd(zeros(1,4),R,N);

%% Initialize variables
sig = zeros(N,1);
x = zeros(size(A{1},1),N);
y = zeros(size(C{1},1),N);
u = ones(1,N);

%% Simulate model
sig = zeros(N,1);
sig(1:round(N/3)) = 1;
sig(round(N/3):round(2*N/3)) = 2;
sig(round(2*N/3):end) = 3;

% sig(1) =  randsrc(1,1,[1,2,3;rho0(1),rho0(2),rho0(3)]); %Random variable with value 1 or 2 with prob of 0.5
x(:,1) = mvnrnd(x0,P0);
% sys_mats.x0 = x(1,:);
for k = 1:N
    if(k<N)
        x(:,k+1) = A{sig(k)}*x(:,k)+B{sig(k)}*u(:,k)+G{sig(k)}*w(k,:)';
        
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
