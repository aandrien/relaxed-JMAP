function [x,y,y_no_noise,sig,sys_mats] = example_blind_seq_est(N,Nf)
%Example 2 from 'An online sequential algorithm for the estimation of
%transition probabilities for jump Markov linear systems 2006 Orguner'

%% Main code
%% Settings
% N = end time
% Nf = 5; %number of filter coefficients

%% System
A = eye(Nf);
B = zeros(Nf,1);
% C = zeros(1,Nf);
C = cell(1,N);
G = eye(Nf);
D = 0;
H = 1;
M = [-1,1];

sys_mats.A = A; sys_mats.B = B; sys_mats.G = G;
sys_mats.D = D; sys_mats.H = H; sys_mats.M = M;

x0 = [0.620; 0.560; 0.480; 0.460; 0.220];
% x0 = rand(Nf,1);
P0 = eye(Nf)*1e-1;
Q = eye(Nf)*1e-6;
R = 1e-0;

rho0(1) = 0.5;
rho0(2) = 0.5;

rho_mat = [0.5 0.5;
    0.5 0.5];

if(any(sum(rho_mat,1)~=1))
    disp('Error, probabilities do not sum up to one!')
    return
end

sys_mats.x0 = x0; sys_mats.P0 = P0; sys_mats.Q = Q; sys_mats.R = R;
sys_mats.rho0 = rho0; sys_mats.rho = rho_mat;

w = mvnrnd(zeros(Nf,1),Q,N);
% w = zeros(N,Nf);
v = mvnrnd(0,R,N);

%% Initialize variables
sig = zeros(N,1);
x = zeros(size(A,1),N);
y = zeros(size(C,1),N);
y_no_noise = zeros(size(C,1),N);

%% Simulate model
sig(1) =  randsrc(1,1,[1,2;rho0(1),rho0(2)]); %Random variable with value 1 or 2 with prob of 0.5
x(:,1) = mvnrnd(x0,P0);
% sys_mats.x0 = x(:,1);
% x(:,1) = x0;
C{1} = ones(1,Nf);
C{1}(1) = M(sig(1));
for k = 1:N
    y(:,k) = C{k}*x(:,k)+H*v(k);
    y_no_noise(:,k) = C{k}*x(:,k);
    
    if(k<N)
%         x(:,k+1) = A*x(:,k)+G*w(k,:)';
        x(:,k+1) = A*x(:,k);
        
        if(sig(k) == 1)
            sig(k+1) = randsrc(1,1,[1,2;rho_mat(1,1),rho_mat(2,1)]); %Random variable 1 with prob 0.6 and 2 with prob 0.4
        elseif(sig(k) == 2)
            sig(k+1) = randsrc(1,1,[1,2;rho_mat(1,2),rho_mat(2,2)]); %Random variable 1 with prob 0.15 and 2 with prob 0.85
        end
        C{k+1} = circshift(C{k},1);
        C{k+1}(1) = M(sig(k+1));
    end   
end
 sys_mats.C = C;
SNR = snr(y,v')
yn = awgn( y_no_noise,SNR,'measured');
end
