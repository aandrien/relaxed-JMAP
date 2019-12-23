function xs = MJLS_smooth_using_xpred2(x,xpred,P,Ppred,sig)

global A B G N Q

xs(:,N) = x(:,N);
for k = N-1:-1:1
    Gj = G{sig(k+1)};
    Aj = A{sig(k+1)};
    Bj = B{sig(k+1)};
%     Qj = Gj*Q*Gj';
    Qj = Q;
    C = P(:,:,k)*Aj'*inv(Ppred(:,:,k+1));
    xs(:,k) = x(:,k)+C*(xs(:,k+1)-xpred(:,k+1));
end
