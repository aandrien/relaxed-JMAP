function [xi,P,xi_pred] = MJLS_KF_BSE_EM(xi_prev,Pprev,sig,yk1,uk,uk1)

% Kalman filter form paper Expectation Maximization Algorithms for MAP 
% Estimation of Jump Markov Linear Systems 1999 Logothetis

global A B C D G H Q R ell

Ck = flipud(sig)';

%% Run filter
nx = size(A,1);
nu = size(B,2);
nw = size(G,2);
ny = size(Ck,1);

Astar = [A, zeros(nx,nx);
    eye(nx), zeros(nx,nx)];

Bstar = [B, zeros(nx,nu);
    zeros(nx,nu), zeros(nx,nu)];

Gstar = [G, zeros(nx,nw);
    zeros(nx,nw), zeros(nx,nw)];

Cstar = [Ck zeros(ny,nx)];

Dstar = [D zeros(ny,nu)];

ustar = [uk1;uk];

Qstar = [Q zeros(nw,nw);
        zeros(nw,nw), Q];

xi_pred = Astar*xi_prev+Bstar*ustar;
Ppred = Astar*Pprev*Astar'+Gstar*Qstar*Gstar';
z = Cstar*xi_pred+Dstar*ustar;
S = Cstar*Ppred*Cstar'+H*R*H';
xi = xi_pred+Ppred*Cstar'*inv(S)*(yk1-z);
P = Ppred-Ppred*Cstar'*inv(S)*Cstar*Ppred;




% Dk = D;
% Hk = H;
% Gk = G;
% Ak = A;
% Bk = B;
% Rk = Hk*R*Hk';
% Qk = Gk*Q*Gk';
% 
% xpred = Ak*xprev+Bk*uk;
% Ppred = Ak*Pprev*Ak'+Gk*Qk*Gk'/(1+e);
% L = inv(Rk/(1+e)+Ck*Ppred*Ck');
% K = Ppred*Ck'*L;
% ypred = Ck*xpred+Dk*uk1;
% yerr = yk1-ypred;
% x = xpred+K*yerr;
% P = Ppred-K*Ck*Ppred;



end