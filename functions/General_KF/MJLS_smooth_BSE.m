function xs = MJLS_smooth_BSE(x,P,u,sig)

global A B G Q

N = size(x,2);
xs(:,N) = x(:,N);

Gj = G;
Aj = A;
Bj = B;
Qj = Gj*Q*Gj';
if(all(Qj == 0))
    invQj = zeros(size(Q,1));
else
    invQj = inv(Qj);
end

for j = N-1:-1:1
    what = -inv(Gj'*inv(Aj)'*inv(P(:,:,j))*inv(Aj)*Gj+invQj)*Gj'*inv(Aj)'*inv(P(:,:,j))*(x(:,j)-inv(Aj)*(xs(:,j+1)-Bj*u(:,j)));
    xs(:,j) = inv(Aj)*(xs(:,j+1)-Gj*what-Bj*u(:,j));
end
