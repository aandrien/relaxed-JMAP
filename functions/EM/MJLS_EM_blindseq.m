function [x_out,sig_out,n_it] = MJLS_EM_blindseq(Nf,y,u,sig_in,stopcrit,x)

%% Global variables
global A B C D G H x0 P0 Q R s N n rho0 rho M

xs_old = zeros(n,N);
% sig_old = randsrc(N,1,[1,2,3;0.3,0.3,0.4])
sig_old = 5*ones(N,1);
xKFs = x;


sig = sig_in;
% sig = randsrc(N,1,[1,-1;0.5,0.5]);
% sig = [2;2;1;2;1];
n_it = 0;
% while(norm(xs_old-xs,2)>stopcrit || norm(sig_old-sig)>stopcrit || n_it < 10)
while((norm(xs_old-xKFs,2)>stopcrit || norm(sig_old-sig)>stopcrit) && n_it<50)
    % while(norm(sig_old-sig)>stopcrit && n_it<100)
    % while(n_it<100)
    %     for j = 1
    xs_old = xKFs;
    sig_old = sig;
    
    %% Initialization
    P_KF = zeros(n,n,N);
    xKF = zeros(n,N);
    
    %% Part 1: Maximization using Kalman Smoother
    % Time 0 KF
    C0 = ones(1,Nf);
    D0 = D;
    H0 = H;
    R0 = H0*R*H0';
    L0 = inv(R0+C0*P0*C0');
    K0 = P0*C0'*L0;
    xKF(:,1) = x0+K0*(y(1)-C0*x0-D0*u(1));
    P_KF(:,:,1) = P0-K0*C0*P0;
    
    % Run KF from k = 1,...,N
    Nt = 1;
    for k = 1:Nf-1
        Nt = Nt+1;
        [xKF(:,k+1),P_KF(:,:,k+1)] = MJLS_KF_BSE(xKF(:,k),P_KF(:,:,k),[ones(Nf-Nt,1);M(sig(k-Nt+2:k+1))'],y(k+1),u(k),u(k+1),0,0);
        
    end
    
    for k = Nf:N-1
        % KF with perfect knowledge of mode
        [xKF(:,k+1),P_KF(:,:,k+1)] = MJLS_KF_BSE(xKF(:,k),P_KF(:,:,k),M(sig(k-Nf+2:k+1))',y(k+1),u(k),u(k+1),0,0);
    end
    
    % Smoothing step
    xKFs = MJLS_smooth_BSE(xKF,P_KF,u,sig);
    
    %% Part 2: Maximization using Viterbi Algorithm
    % Calculate b_i(y_k,x_k^(l))
    b = zeros(s,N);
    Ri = H*R*H';
    for ii = 1:s
        %         Ci = [M(ii),ones(1,Nf-1)];
%                 Ci = repmat(M(ii),1,5);
        Ci = ones(1,Nf);
        Ci(1) = M(ii);
        %         b(ii,1) = -0.5*mnorm(y(:,1)-Ci*xKFs(:,1),inv(Ri))^2;
        for k = 1:N
            b(ii,k) = -0.5*mnorm(y(:,k)-Ci*xKFs(:,k),inv(Ri))^2;
%             b(ii,k) = -0.5*(y(:,k)-Ci*xKFs(:,k))'*inv(Ri)*(y(:,k)-Ci*xKFs(:,k));
            Ci = circshift(Ci,1);
            if(k<N)
                Ci(1) = M(ii);
            end
        end
    end
    
    % Initialization
    delta = zeros(s,N);
    psi_t = zeros(s,N);
    
    delta_ = zeros(s,N);
    psi_t_ = zeros(s,N);
    
    for ii = 1:s
        delta(ii,1) = log(rho0(ii))+b(ii,1);
        delta_(ii,1) = log(rho0(ii))+b(ii,1);
    end
    
    % Recursion
    for k = 2:N
        for jj = 1:s
            delta(jj,k) = b(jj,k)+max(delta(:,k-1)'+log(rho(:,jj))');
            [~,psi_t(jj,k)] = max(delta(:,k-1)'+log(rho(:,jj))');
        end
    end
    
    for k = 2:N
        for jj = 1:s
            delta_(jj,k) = max(b(:,k)'+delta_(:,k-1)'+log(rho(:,jj))');
            [~,psi_t_(jj,k)] = max(b(:,k)'+delta_(:,k-1)'+log(rho(:,jj))');
        end
    end
    
    % Termination
    sig = zeros(N,1);
    [~,sig(N,1)] = max(delta(:,N));
    
    % Backtracking
    for k = N-1:-1:1
        sig(k,1) = psi_t(sig(k+1,1),k+1);
    end
    n_it = n_it+1;
    
end
x_out = xKFs;
sig_out = sig;
end