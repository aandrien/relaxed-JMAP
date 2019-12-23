function xs = MJLS_smooth_parfor(sys,x,P,u,sig)

A = sys.A; B = sys.B; G = sys.G;
C = sys.C; D = sys.D; H = sys.H;
x0 = sys.x0; P0 = sys.P0; Q = sys.Q; R = sys.R;
rho0 = sys.rho0; rho = sys.rho; N = sys.N;

ell0 = -log(rho0);
ell = -log(rho);

n = size(A{1},1);
s = size(A,2);


xs(:,N) = x(:,N);
for j = N-1:-1:1
    Gj = G{sig(j)};
    Aj = A{sig(j)};
    Bj = B{sig(j)};
    Qj = Gj*Q*Gj';
    what = -inv(Gj'*inv(Aj)'*inv(P(:,:,j))*inv(Aj)*Gj+inv(Qj))*Gj'*inv(Aj)'*inv(P(:,:,j))*(x(:,j)-inv(Aj)*(xs(:,j+1)-Bj*u(:,j)));
    xs(:,j) = inv(Aj)*(xs(:,j+1)-Gj*what-Bj*u(:,j));
end
