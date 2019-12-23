function [x_RDP,P_RDP,cRDP] = RDP_post(sys,Pset,y,u)

n = sys.n;
N = sys.N;

x_RDP = zeros(n,N+1);
P_RDP = zeros(n,n,N+1);
cRDP = zeros(N+1,1);

ii = Pset(1,1);
x_RDP(:,1) = sys.x0;
P_RDP(:,:,1) = sys.P0;
cRDP(1,1) = sys.ell0(ii);

for k = 1:N
    % KF with perfect knowledge of mode
    [x_RDP(:,k+1),P_RDP(:,:,k+1),cRDP(k+1)] = MJLS_KF2(sys,x_RDP(:,k),P_RDP(:,:,k),Pset(k,1),Pset(k+1,1),y(:,k+1),u(:,k),u(:,k+1),cRDP(k),0);
end
