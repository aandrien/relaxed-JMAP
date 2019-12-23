function xis = MJLS_smooth_BSE_EM(xi,xi_pred,P,u,y,sig,Nf)

global A B D G H Q R C

N = size(xi,2);
xis(:,N) = xi(:,N);



%% Run filter
nx = size(A,1);
nu = size(B,2);
nw = size(G,2);
ny = size(1,1);

Astar = [A, zeros(nx,nx);
    eye(nx), zeros(nx,nx)];

Ck = flipud(sig(end-Nf+1:end))';
Cstar = [Ck zeros(ny,nx)];



lambda(:,N) = zeros(2*nx,1);
J(:,:,N) = Cstar'*inv(H*R*H')*Cstar;
Nt = Nf;
for k = N:-1:1
    xis(:,k) = xi(:,k)-P(:,:,k)*Astar'*lambda(:,k);
    
    if(k>Nf)
        Ck = flipud(sig(k-Nf+1:k))';
    else
        Ck = flipud([ones(Nf-Nt,1);sig(k-Nt+1:k)])';
        Nt = Nt-1;
    end
    Cstar = [C{k} zeros(ny,nx)];
    
    J(:,:,k) = Cstar'*inv(H*R*H')*Cstar;
    if(k>1)
        lambda(:,k-1) = (eye(2*nx)-P(:,:,k)*J(:,:,k))'* ...
            (Astar'*lambda(:,k)-Cstar'*inv(H*R*H')*(y(:,k)-Cstar*xi_pred(:,k)));
    end
end
