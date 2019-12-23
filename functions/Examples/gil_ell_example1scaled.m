function [x,y,sig,sys_mats] = gil_ell_example1scaled(N,T)
%Example 1 from 'An online sequential algorithm for the estimation of
%transition probabilities for jump Markov linear systems 2006 Orguner'

%% Main code
%% Settings
% N = end time

%% System
A{1} = 1; A{2} = 1;
B{1} = 0; B{2} = 0;
G{1} = 1; G{2} = 1;
C{1} = 1; C{2} = 1;
D{1} = 0; D{2} = 0;
H{1} = 20; H{2} = 1;
% H{1} = 10; H{2} = 10;

% Rescale
A{1} = inv(T)*A{1}*T ; A{2} = inv(T)*A{2}*T;
B{1} = inv(T)*B{1}; B{2} = inv(T)*B{2};
G{1} = inv(T)*G{1}; G{2} = inv(T)*G{2};

C{1} = inv(T)*C{1}*T; C{2} = inv(T)*C{2}*T;
D{1} = inv(T)*D{1}; D{2} = inv(T)*D{2};
H{1} = inv(T)*H{1}; H{2} = inv(T)*H{2};


sys_mats.A = A; sys_mats.B = B; sys_mats.G = G;
sys_mats.C = C; sys_mats.D = D; sys_mats.H = H;

x0 = 0;
P0 = (inv(T)*2)^2;
Q = 10^2;
R = 1^2;

rho0(1) = 0.5;
rho0(2) = 0.5;

% rho = [0.6, 0.85
%     0.4, 0.15];
rho = [0.15, 0.4
    0.85, 0.6];


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

%% Generate noise
w = mvnrnd(0,Q,N);
v = mvnrnd(0,R,N);

%% Initialize variables
sig = zeros(N,1);
x = zeros(size(A,1),N);
y = zeros(size(C,1),N);

%% Simulate model
sig(1) =  randsrc(1,1,[1,2;rho0(1),rho0(2)]); %Random variable with value 1 or 2 with prob of 0.5
x(1) = mvnrnd(x0,P0);
% sys_mats.x0 = x(1);

for k = 2:N+1
    if(sig(k-1) == 1)
        sig(k) = randsrc(1,1,[1,2;rho(1,1),rho(2,1)]); %Random variable 1 with prob 0.6 and 2 with prob 0.4
    elseif(sig(k-1) == 2)
        sig(k) = randsrc(1,1,[1,2;rho(1,2),rho(2,2)]); %Random variable 1 with prob 0.15 and 2 with prob 0.85
    end
    x(k) = A{sig(k)}*x(k-1)+G{sig(k)}*w(k-1);
    y(k) = C{sig(k)}*x(k)+H{sig(k)}*v(k-1);
end