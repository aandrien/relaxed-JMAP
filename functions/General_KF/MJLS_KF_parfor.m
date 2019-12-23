function [x,P,c] = MJLS_KF_parfor(sys,xprev,Pprev,sigk,sigk1,yk1,uk,uk1,cprev,e)

if(nargin<10)
    e = 0;
end

A = sys.A; B = sys.B; G = sys.G;
C = sys.C; D = sys.D; H = sys.H;
Q = sys.Q; R = sys.R;
rho0 = sys.rho0; rho = sys.rho;

ell = -log(rho);

%% Run filter
Ck = C{sigk1};
Dk = D{sigk1};
Hk = H{sigk1};
Gk = G{sigk1};
Ak = A{sigk1};
Bk = B{sigk1};
Rk = Hk*R*Hk';
Qk = Gk*Q*Gk';

xpred = Ak*xprev+Bk*uk;
Ppred = Ak*Pprev*Ak'+Qk/(1+e);
L = inv(Rk/(1+e)+Ck*Ppred*Ck');
K = Ppred*Ck'*L;
ypred = Ck*xpred+Dk*uk1;
yerr = yk1-ypred;
x = xpred+K*yerr;
P = Ppred-K*Ck*Ppred;
c = 0.5*mnorm(yerr,L)^2+(1+e)*ell(sigk1,sigk)+cprev;

end