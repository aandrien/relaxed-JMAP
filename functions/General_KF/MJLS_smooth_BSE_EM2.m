function xis = MJLS_smooth_BSE_EM2(xi,xi_pred,P,Ppred)

global A

N = size(xi,2);
xis(:,N) = xi(:,N);

%% Run filter
nx = size(A,1);

Astar = [A, zeros(nx,nx);
    eye(nx), zeros(nx,nx)];

for k = N-1:-1:1
    CCk = P(:,:,k)*Astar'*inv(Ppred(:,:,k));
    xis(:,k) = xi(:,k)+CCk*(xis(:,k+1)-xi_pred(:,k+1));
end
