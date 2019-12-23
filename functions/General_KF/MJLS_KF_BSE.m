function [x,P,c] = MJLS_KF_BSE(xprev,Pprev,sig,yk1,uk,uk1,cprev,e)

global A B C D G H Q R ell

if(nargin<9)
    e = 0;
end

Ck = flipud(sig)';

%% Run filter
Dk = D;
Hk = H;
Gk = G;
Ak = A;
Bk = B;
Rk = Hk*R*Hk';
Qk = Gk*Q*Gk';

xpred = Ak*xprev+Bk*uk;
Ppred = Ak*Pprev*Ak'+Gk*Qk*Gk'/(1+e);
L = inv(Rk/(1+e)+Ck*Ppred*Ck');
K = Ppred*Ck'*L;
ypred = Ck*xpred+Dk*uk1;
yerr = yk1-ypred;
x = xpred+K*yerr;
P = Ppred-K*Ck*Ppred;
c = 0.5*mnorm(yerr,L)^2+(1+e)*ell(1,1)+cprev;


end