function J = cost_checktmp(xs,u,y,sig,e)

global N A B C D G H Q R x0 P0 ell0 ell

J = zeros(N,1);

H0 = H{sig(1)};
R0 = H0*R*H0';
J(1) = 0.5*(mnorm(xs(:,1)-x0,inv(P0))^2+mnorm(y(:,1)-C{sig(1)}*xs(:,1)-D{sig(1)}*u(1),inv(R0))^2)+ell0(sig(1));

for kk = 2:N
    if(size(G{sig(kk-1)},1)==size(G{sig(kk-1)},2))
        wk = inv(G{sig(kk-1)})*(xs(:,kk)-A{sig(kk-1)}*xs(:,kk-1)-B{sig(kk-1)}*u(kk-1));
    else
        wk = pinv(G{sig(kk-1)})*(xs(:,kk)-A{sig(kk-1)}*xs(:,kk-1)-B{sig(kk-1)}*u(kk-1));
    end
    Hk = H{sig(kk)};
    Gk = G{sig(kk-1)};
    Rk = Hk*R*Hk';
    Qk = Q;
    J(kk) = 0.5*(mnorm(wk,inv(Qk))^2+mnorm(y(:,kk)-C{sig(kk)}*xs(:,kk)-D{sig(kk)}*u(kk),inv(Rk))^2)+ell(sig(kk),sig(kk-1));
    if(kk == N)
        J(kk) = (1+e)*(0.5*(mnorm(wk,inv(Qk))^2+mnorm(y(:,kk)-C{sig(kk)}*xs(:,kk)-D{sig(kk)}*u(kk),inv(Rk))^2)+ell(sig(kk),sig(kk-1)));
    end
end