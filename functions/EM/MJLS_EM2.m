function [x_out,sig_out,n_it] = MJLS_EM(Nf,y,u,sig_in,stopcrit,x)

%% Global variables
global A B C D G H x0 P0 Q R s N n rho0 rho M

xs_old = zeros(n,N);
% sig_old = randsrc(N,1,[1,2,3;0.3,0.3,0.4])
sig_old = 5*ones(N,1);
xs = x;

sig = sig_in;
% sig = randsrc(N,1,[1,-1;0.5,0.5]);
% sig = [2;2;1;2;1];
n_it = 0;
% while(norm(xs_old-xs,2)>stopcrit || norm(sig_old-sig)>stopcrit || n_it < 10)
% while(norm(xs_old-xs,2)>stopcrit || norm(sig_old-sig)>stopcrit)
while(norm(sig_old-sig)>stopcrit && n_it<100)
% while(n_it<100)
%     for j = 1
    xs_old = xs;
    sig_old = sig;
    
    %% Initialization
    %     P = zeros(n,n,N);
    %     P(:,:,1) = eye(n);
    %     x = zeros(n,N);
    
    %% Part 1: Maximization using Kalman Smoother
    % Run KF from k = 1,...,N
    xi(:,1) = [x0;x0];
    xi_pred(:,1) = xi(:,1);
    % P(:,:,1) = [P0 P0;P0 P0];
    P(:,:,1) = blkdiag(P0,P0);
    Ppred(:,:,1) = P(:,:,1);
    Nt = 1;
    for k = 1:Nf-1
        [xi(:,k+1),P(:,:,k+1),xi_pred(:,k+1),Ppred(:,:,k+1)] = MJLS_KF_BSE_EM(xi(:,k),P(:,:,k),[ones(Nf-Nt,1);sig(k-Nt+1:k)],y(k),u(k),u(k+1));
        Nt = Nt+1;
    end
    
    for k = Nf:N-1
        [xi(:,k+1),P(:,:,k+1),xi_pred(:,k+1),Ppred(:,:,k+1)] = MJLS_KF_BSE_EM(xi(:,k),P(:,:,k),sig(k-Nf+1:k),y(k),u(k),u(k+1));
    end
    
    % Smoothing step
    xis = MJLS_smooth_BSE_EM2(xi,xi_pred,P,Ppred);
    xs = xis(1:Nf,:);
    
    %% Part 2: Maximization using Viterbi Algorithm
    % Calculate b_i(y_k,x_k^(l))
    b = zeros(s,N);
    Ri = H*R*H';
    for ii = 1:s
%         Ci = [M(ii),ones(1,Nf-1)];
        Ci = repmat(M(ii),1,5);

        b(ii,1) = -0.5*mnorm(y(:,1)-Ci*xs(:,1),inv(Ri))^2;
        for k = 2:N
            b(ii,k) = -0.5*mnorm(y(:,k)-Ci*xs(:,k),inv(Ri))^2;
        end
    end
    
    % Initialization
    delta = zeros(s,N);
    psi_t = zeros(s,N);
    
    for ii = 1:s
        delta(ii,1) = log(rho0(ii))+b(ii,1);
    end
    
    % Recursion
    for k = 2:N
        for jj = 1:s
            delta(jj,k) = b(jj,k)+max(delta(:,k-1)'+log(rho(:,jj))');
            [~,psi_t(jj,k)] = max(delta(:,k-1)'+log(rho(:,jj))');
        end
    end
    
    % Termination
    sig = zeros(N,1);
    [~,sig(N,1)] = max(delta(:,N));
    
    % Backtracking
    for k = N-1:-1:1
        sig(k,1) = psi_t(sig(k+1,1),k+1);
    end
    Mnew = [1,-1];
    sig = M(sig)';
    n_it = n_it+1;
    
end
x_out = xs;
sig_out = sig;
end