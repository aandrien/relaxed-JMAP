function J = cost_check2(sys,xs,u,y,sig)

%% Load Variables
N = sys.N;
A = sys.A;
B = sys.B;
G = sys.G;
C = sys.C;
D = sys.D;
H = sys.H;
R = sys.R;
Q = sys.Q;
x0 = sys.x0;
P0 = sys.P0;
ell0 = sys.ell0;
ell = sys.ell;
n = sys.n;
r = sys.r;

J = zeros(N,1);

%% Calculate Cost
J(1) = 0.5*(mnorm(xs(:,1)-x0,inv(P0)))^2+ell0(sig(1));

for kk = 2:N+1
    Gk = G{sig(kk)};
    if(size(Gk,1)==size(Gk,2))
        wk = (xs(:,kk)-A{sig(kk)}*xs(:,kk-1)-B{sig(kk)}*u(kk-1));
    else
        wk = (xs(:,kk)-A{sig(kk)}*xs(:,kk-1)-B{sig(kk)}*u(kk-1));
    end
    Hk = H{sig(kk)};
    Rk = Hk*R*Hk';
    Qk = Gk*Q*Gk';
    J(kk) = 0.5*(mnorm(wk,inv(Qk))^2 ...
        +mnorm(y(:,kk)-C{sig(kk)}*xs(:,kk)-D{sig(kk)}*u(kk),inv(Rk))^2) ...
        +ell(sig(kk),sig(kk-1)) ...
        +(log(2*pi)*n/2+log(det(Qk))/2) ...
        +(log(2*pi)*r/2+log(det(Rk))/2);
end