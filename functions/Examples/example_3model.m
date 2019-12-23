function [x,y,sig,sys_mats] = example_3model(N)
%Example 2 from 'An online sequential algorithm for the estimation of
%transition probabilities for jump Markov linear systems 2006 Orguner'

%% Main code
%% Settings
% N = end time

%% System
A{1} = 0.8; A{2} = 0.9; A{3} = 1;
B{1} = 0; B{2} = 0; B{3} = 0;
G{1} = 1; G{2} = 1; G{3} = 1;
C{1} = 1; C{2} = 2; C{3} = 4;
D{1} = 0; D{2} = 0; D{3} = 0;
H{1} = 1; H{2} = 1; H{3} = 1;
% H{1} = 10; H{2} = 10;
sys_mats.A = A; sys_mats.B = B; sys_mats.G = G;
sys_mats.C = C; sys_mats.D = D; sys_mats.H = H;

x0 = 10;
P0 = 2^2;
Q = 2^2;
R = 1;

rho0(1) = 0.2;
rho0(2) = 0.6;
rho0(3) = 0.2;

rho_mat = [0.2 0.25 0.1;
    0.4 0.5 0.1;
    0.4 0.25 0.8];

if(any(sum(rho_mat,1)~=1))
    disp('Error, probabilities do not sum up to one!')
    return
end

sys_mats.x0 = x0; sys_mats.P0 = P0; sys_mats.Q = Q; sys_mats.R = R;
sys_mats.rho0 = rho0; sys_mats.rho = rho_mat;

w = mvnrnd(0,Q,N);
v = mvnrnd(0,R,N);

%% Initialize variables
sig = zeros(N,1);
x = zeros(size(A,1),N);
y = zeros(size(C,1),N);

%% Simulate model
sig(1) =  randsrc(1,1,[1,2,3;rho0(1),rho0(2),rho0(3)]); %Random variable with value 1 or 2 with prob of 0.5
x(1) = mvnrnd(x0,P0);
sys_mats.x0 = x(1);
for k = 1:N
    if(k<N)
        x(k+1) = A{sig(k)}*x(k)+G{sig(k)}*w(k);
        
        if(sig(k) == 1)
            sig(k+1) = randsrc(1,1,[1,2,3;rho_mat(1,1),rho_mat(2,1),rho_mat(3,1)]); %Random variable 1 with prob 0.6 and 2 with prob 0.4
        elseif(sig(k) == 2)
            sig(k+1) = randsrc(1,1,[1,2,3;rho_mat(1,2),rho_mat(2,2),rho_mat(3,2)]); %Random variable 1 with prob 0.15 and 2 with prob 0.85
        elseif(sig(k) == 3)
            sig(k+1) = randsrc(1,1,[1,2,3;rho_mat(1,3),rho_mat(2,3),rho_mat(3,3)]); %Random variable 1 with prob 0.15 and 2 with prob 0.85
        end
    end
    y(k) = C{sig(k)}*x(k)+H{sig(k)}*v(k);
end
