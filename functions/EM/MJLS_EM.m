function [x_out,sig_out,n_it] = MJLS_EM(sys,y,u,sig_in,stopcrit,x)

n = sys.n;
s = sys.s;
N = sys.N;
rho = sys.rho;
rho0 = sys.rho0;
x0 = sys.x0;
P0 = sys.P0;

%% Global variables
xs_old = 100*ones(n,N+1);
% sig_old = randsrc(N,1,[1,2,3;0.3,0.3,0.4])
sig_old = 5*ones(N+1,1);
% xKFs = x;
xKFs = zeros(n,N+1);

% sig = sig_in;
% sig(randi(N+1,3,1))= randi(s,3,1);

% sig = randsrc(N+1,1,[1,2,3;0.4,0.4,0.2]);
sig = randi(s,N+1,1);
% sig = [2;2;1;2;1];
siginit = sig;
n_it = 0;
while((norm(xs_old-xKFs,2)>stopcrit || norm(sig_old-sig)>stopcrit) && n_it<50)
    xs_old = xKFs;
    sig_old = sig;

    %% Initialization
    P_KF = zeros(n,n,N+1);
    xKF = zeros(n,N+1);
    
    %% Part 1: Maximization using Kalman Smoother
    % Time 0 KF
    xKF(:,1) = x0;
    P_KF(:,:,1) = P0;
    
    % Run KF from k = 1,...,N  
    for k = 1:N
        % KF with perfect knowledge of mode
        [xKF(:,k+1),P_KF(:,:,k+1),~] = MJLS_KF2(sys,xKF(:,k),P_KF(:,:,k),sig(k),sig(k+1),y(:,k+1),u(k),u(k+1),0);
    end

    % Smoothing step
    xKFs = MJLS_smooth2(sys,xKF,P_KF,u,sig);
    J(:,n_it+1) = cost_check2(sys,xKFs,u,y,sig);
    %% Part 2: Maximization using Viterbi Algorithm
    % Calculate b_i(y_k,x_k^(l))
    b = zeros(s,N+1);
    
    for ii = 1:s
        Ci = sys.C{ii};
        Ri = sys.H{ii}*sys.R*sys.H{ii}';
        Qi = sys.G{ii}*sys.Q*sys.G{ii}';
        Ai = sys.A{ii};
        b(ii,1) = -0.5*mnorm(xKFs(:,1)-x0,inv(P0))^2;
        for k = 2:N+1
%             b(ii,k) = -0.5*mnorm(y(:,k)-Ci*xKFs(:,k),inv(Ri))^2;
            b(ii,k) = -0.5*mnorm(y(:,k)-Ci*xKFs(:,k),inv(Ri))^2 ...
            -0.5*mnorm(xKFs(:,k)-Ai*xKFs(:,k-1),inv(Qi))^2; ...
%             -log(det(Qi))/2 ...
%             -log(det(Ri))/2;
        end
    end
    
    % Initialization
    delta = zeros(s,N+1);
    psi_t = zeros(s,N+1);
    
    for ii = 1:s
        delta(ii,1) = log(rho0(ii))+b(ii,1);
    end
    
    % Recursion
    for k = 2:N+1
        for jj = 1:s
            delta(jj,k) = b(jj,k)+max(delta(:,k-1)'+log(rho(:,jj))');
            [~,psi_t(jj,k)] = max(delta(:,k-1)'+log(rho(:,jj))');
        end
    end
    
    % Termination
    sig = zeros(N+1,1);
    [~,sig(N+1,1)] = max(delta(:,N+1));
    
    % Backtracking
    for k = N:-1:1
        sig(k,1) = psi_t(sig(k+1,1),k+1);
    end
    n_it = n_it+1;

end
x_out = xKFs;
sig_out = sig;
end