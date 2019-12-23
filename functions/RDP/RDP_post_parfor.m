function [x_RDP,P_RDP,cRDP] = RDP_post_parfor(sys,Pset,y,u)

A = sys.A; B = sys.B; G = sys.G;
C = sys.C; D = sys.D; H = sys.H;
x0 = sys.x0; P0 = sys.P0; Q = sys.Q; R = sys.R;
rho0 = sys.rho0; rho = sys.rho; N = sys.N;

ell0 = -log(rho0);
ell = -log(rho);

n = size(A{1},1);
s = size(A,2);

x_RDP = zeros(n,N);
P_RDP = zeros(n,n,N);
cRDP = zeros(N,1);

ii = Pset(1,1);
C0 = C{ii};
D0 = D{ii};
H0 = H{ii};
R0 = H0*R*H0';
L0 = inv(R0+C0*P0*C0');
K0 = P0*C0'*L0;
x_RDP(:,1) = x0+K0*(y(:,1)-C0*x0-D0*u(:,1));
P_RDP(:,:,1) = P0-K0*C0*P0;
cRDP(1,1) = 0.5*mnorm(y(:,1)-C0*x0-D0*u(:,1),L0)^2+ell0(ii);

for k = 1:N-1
    % KF with perfect knowledge of mode
    [x_RDP(:,k+1),P_RDP(:,:,k+1),cRDP(k+1)] = MJLS_KF_parfor(sys,x_RDP(:,k),P_RDP(:,:,k),Pset(k,1),Pset(k+1,1),y(:,k+1),u(:,k),u(:,k+1),cRDP(k),0);
end
