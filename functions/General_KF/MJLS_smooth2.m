function xs = MJLS_smooth2(sys,x,P,u,sig)

N = sys.N;

xs(:,N+1) = x(:,N+1);
for j = N:-1:1
    Gj = sys.G{sig(j+1)};
    Aj = sys.A{sig(j+1)};
    Bj = sys.B{sig(j+1)};
%     Qj = Gj*sys.Q*Gj';
    Qj = sys.Q;
    what = -inv(Gj'*inv(Aj)'*inv(P(:,:,j))*inv(Aj)*Gj+inv(Qj))*Gj'*inv(Aj)'*inv(P(:,:,j))*(x(:,j)-inv(Aj)*(xs(:,j+1)-Bj*u(:,j)));
    xs(:,j) = inv(Aj)*(xs(:,j+1)-Gj*what-Bj*u(:,j));
end
