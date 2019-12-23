function J = cost_check_parfor(sys,xs,u,y,sig)

A = sys.A; B = sys.B; G = sys.G;
C = sys.C; D = sys.D; H = sys.H;
x0 = sys.x0; P0 = sys.P0; Q = sys.Q; R = sys.R;
rho0 = sys.rho0; rho = sys.rho; N = sys.N;

ell0 = -log(rho0);
ell = -log(rho);

n = size(A{1},1);
s = size(A,2);


J = zeros(N,1);

H0 = H{sig(1)};
R0 = H0*R*H0';
J(1) = 0.5*(mnorm(xs(:,1)-x0,inv(P0))^2+mnorm(y(1)-C{sig(1)}*xs(:,1)-D{sig(1)}*u(1),inv(R0))^2)+ell0(sig(1));

for kk = 2:N
    wk = inv(G{sig(kk-1)})*(xs(:,kk)-A{sig(kk-1)}*xs(:,kk-1)-B{sig(kk-1)}*u(kk-1));
    Hk = H{sig(kk)};
    Gk = G{sig(kk-1)};
    Rk = Hk*R*Hk';
    Qk = Gk*Q*Gk';
    J(kk) = 0.5*(mnorm(wk,inv(Qk))^2+mnorm(y(kk)-C{sig(kk)}*xs(:,kk)-D{sig(kk)}*u(kk),inv(Rk))^2)+ell(sig(kk),sig(kk-1));
end