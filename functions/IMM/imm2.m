function [xhat,Phat,x,P,mu_] = imm2(sys,xprev,Pprev,mu_p,yk1,uk,uk1)

s = sys.s;
n= sys.n;

%% Calculation of mixing probabilities
mu_ = zeros(s,s);
p = sys.rho';
cbar = p'*mu_p;
for ii = 1:s
    for jj = 1:s
        mu_(ii,jj) = p(ii,jj)*mu_p(ii)/cbar(jj);
    end
end

%% Mixing initial condition for the filter matched to j
P0j = zeros(n,n,s);
x0 = zeros(n,s);
for jj = 1:s
    % 3) Compute xhat(mod)
    x0(:,jj) = xprev*mu_(:,jj);
    % 4) Compute Sigma(mod).
    
    for ii = 1:s
        xd = xprev(:,ii)-x0(:,jj);
        P0j(:,:,jj) = P0j(:,:,jj)+mu_(ii,jj)*(Pprev(:,:,ii)+xd*xd');
    end
end

%% Perform 1 time step for each of the 's' Kalman filters
x = zeros(n,s);
P = zeros(n,n,s);
Lambda = zeros(s,1);
for jj = 1:s
    [x(:,jj),P(:,:,jj),ypred,L] = MJLS_IMM2(sys,x0(:,jj),P0j(:,:,jj),jj,jj,yk1,uk,uk1,0);
    Lambda(jj) = max(1e-10,exp(-0.5*(yk1-ypred)'*L*(yk1-ypred))/sqrt(2*pi*det(L)));
end

%% Mode probability Update
c = Lambda'*cbar;
mu_ = Lambda.*cbar./c;

%% Output estimate calculation
xhat = x*mu_;
Phat = zeros(n,n);
for jj = 1:s
    Phat = Phat+mu_(jj)*(P(:,:,jj)+(x(:,jj)-xhat)*(x(:,jj)-xhat)');
end

end