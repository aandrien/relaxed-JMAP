function [x,y,u,sig,sys_mats] = example_coloured_noise(N)
%Example of moving target with constant velocity and constant turn rates

%% Main code
%% Settings
% N = end time
c1 = -3.83; c2 = 2.62; c3 = -1.15;
q = 0.9;

%% System
nModes = 2;
% Rest of system is same for each sigma
for jj = 1:nModes
    A{jj} = [0 0 0 0; 
        1 0 0 0;
        0 1 0 0;
        0 0 1 0];
    B{jj} = zeros(4,1);
    G{jj} = [1 0 0 0]'; 
    C{jj} = [1 c1 c2 c3];
    H{jj} = eye(1);
end
D{1} = -q;
D{2} = q;

sys_mats.A = A; sys_mats.B = B; sys_mats.G = G;
sys_mats.C = C; sys_mats.D = D; sys_mats.H = H;

x0 = zeros(4,1);
P0 = 0.01*eye(4);
Q = 0.1.*eye(1);
R = 0.001.*eye(1);

rho0 = [0.9; 0.1];

rho_mat = [0.9 0.1;
    0.1 0.9];

if(any(sum(rho_mat,1)~=1))
    disp('Error, probabilities do not sum up to one!')
    return
end

sys_mats.x0 = x0; sys_mats.P0 = P0; sys_mats.Q = Q; sys_mats.R = R;
sys_mats.rho0 = rho0; sys_mats.rho = rho_mat;

w = mvnrnd(zeros(1,1),Q,N);
v = mvnrnd(zeros(1,1),R,N);

%% Initialize variables
sig = zeros(N,1);
x = zeros(size(A{1},1),N);
y = zeros(size(C{1},1),N);
u = ones(1,N);

%% Simulate model
sig(1) =  randsrc(1,1,[1,2;rho0(1),rho0(2)]); %Random variable with value 1 or 2 with prob of 0.5
x(:,1) = mvnrnd(x0,P0);
% sys_mats.x0 = x(1,:);
for k = 1:N
    if(k<N)
        x(:,k+1) = A{sig(k)}*x(:,k)+B{sig(k)}*u(:,k)+G{sig(k)}*w(k,:)';
        %         x(:,k+1) = A{sig(k)}*x(:,k)+B{sig(k)}*u(:,k);
        
        if(sig(k) == 1)
            sig(k+1) = randsrc(1,1,[1,2;rho_mat(1,1),rho_mat(2,1)]); %Random variable 1 with prob 0.6 and 2 with prob 0.4
        elseif(sig(k) == 2)
            sig(k+1) = randsrc(1,1,[1,2;rho_mat(1,2),rho_mat(2,2)]); %Random variable 1 with prob 0.15 and 2 with prob 0.85
        end
    end
    y(:,k) = C{sig(k)}*x(:,k)+D{sig(k)}*u(:,k)+H{sig(k)}*v(k,:)';
end
