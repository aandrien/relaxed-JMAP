function xs = MJLS_smoothbar(x,P,u,sig,e)

global A B G N Q

xs(:,N) = x(:,N);
for j = N-1:-1:1
    Gj = G{sig(j)};
    Aj = A{sig(j)};
    Bj = B{sig(j)};
    Qj = Gj*Q*Gj';
    if(j==N-1)
        Qj = Q/(1+e);
    end
    what = -inv(Gj'*inv(Aj)'*inv(P(:,:,j))*inv(Aj)*Gj+inv(Qj))*Gj'*inv(Aj)'*inv(P(:,:,j))*(x(:,j)-inv(Aj)*(xs(:,j+1)-Bj*u(:,j)));
    xs(:,j) = inv(Aj)*(xs(:,j+1)-Gj*what-Bj*u(:,j));
end
