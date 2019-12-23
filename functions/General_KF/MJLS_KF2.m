function [x,P,c] = MJLS_KF2(sys,xprev,Pprev,sigk,sigk1,yk1,uk,uk1,cprev,e)

if(nargin<10)
    e = 0;
end

n = sys.n;
r = sys.r;

%% Read Variables
Ck = sys.C{sigk1};
Dk = sys.D{sigk1};
Hk = sys.H{sigk1};
Gk = sys.G{sigk1};
Ak = sys.A{sigk1};
Bk = sys.B{sigk1};
Rk = Hk*sys.R*Hk';
Qk = Gk*sys.Q*Gk';
ell = sys.ell;

%% Run filter
xpred = Ak*xprev+Bk*uk;
Ppred = Ak*Pprev*Ak'+Qk/(1+e);
L = inv(Rk/(1+e)+Ck*Ppred*Ck');
K = Ppred*Ck'*L;
ypred = Ck*xpred+Dk*uk1;
yerr = yk1-ypred;
x = xpred+K*yerr;
P = Ppred-K*Ck*Ppred;
c = 0.5*mnorm(yerr,L)^2 ...
    +(1+e)*ell(sigk1,sigk) ...
    +cprev ...
    +(1+e)*(log(2*pi)*n+log(det(Qk)))/2 ...
    +(1+e)*(log(2*pi)*r+log(det(Rk)))/2;
%     -(1+e)*(log(det(Qk))/2+log(det(Rk))/2);
%     +(1+e)*(log(det(Qk))/2) ...
%     +(1+e)*(log(det(Rk))/2);



end