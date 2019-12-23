function [xout,Pout,cout,ind_min,Pset_out] = MJLS_RDP_LMI_parfor(sys,xprev,Pprev,cprev,yk1,uk,uk1,Pset_prev,e)

A = sys.A; B = sys.B; G = sys.G;
C = sys.C; D = sys.D; H = sys.H;
x0 = sys.x0; P0 = sys.P0; Q = sys.Q; R = sys.R;
rho0 = sys.rho0; rho = sys.rho;

ell0 = -log(rho0);
ell = -log(rho);

n = size(A{1},1);
s = size(A,2);

%% Calculate KF for each entry in set C
cnt = 1;
hst = zeros(size(Pset_prev,2)*s,1);
nxt = zeros(size(Pset_prev,2)*s,1);
x = zeros(n,size(Pset_prev,2)*s);
P = zeros(n,n,size(Pset_prev,2)*s);
c = zeros(size(Pset_prev,2)*s,1);

for ii = 1:size(Pset_prev,2) %Loop over entries in P_{k-1}
    for jj = 1:s %For each entry in P_{k-1} there are s new options (set C_K)
        [x(:,cnt),P(:,:,cnt),c(cnt)] = MJLS_KF_parfor(sys,xprev(:,ii),Pprev(:,:,ii),Pset_prev(end,ii),jj,yk1,uk,uk1,cprev(ii),0);
        hst(cnt) = ii;
        nxt(cnt) = jj;
        cnt = cnt+1;
    end
end

%% RDP procedure using minmax
%[1] initialize set Pset as empty
%[2] Add the element in Cset with the smallest c to Pset
c_hst_nxt = sortrows([c,hst,nxt,[1:cnt-1]'],1); %sort
c = c_hst_nxt(:,1);
hst = c_hst_nxt(:,2);
nxt = c_hst_nxt(:,3);
cnts = c_hst_nxt(:,4);

tokeep = 1;
xopt = x(:,cnts(tokeep));
Pset_tmp = [Pset_prev(:,hst(1));nxt(1)];

%[3] Take one element i in Ck and check if it satisfies relaxed cost
opts = sdpsettings('verbose',0,'solver','sdpt3');

for kk = 2:cnt-1
    %Calculate relaxed cost
    [xbar,Pbar,cbar] = MJLS_KF_parfor(sys,xprev(:,hst(kk)),Pprev(:,:,hst(kk)),Pset_prev(end,hst(kk)),nxt(kk),yk1,uk,uk1,cprev(hst(kk)),e);
    
    %Check condition using LMIs
    nP = size(tokeep,2);
    tau = sdpvar(nP,1);
    PiBar = inv(Pbar);
    barMatrix = [PiBar, -PiBar*xbar; -xbar'*PiBar, xbar'*PiBar*xbar+cbar];
    hatMatrix = zeros(2*n,2*n);
    for jj = 1:nP
        PiHat = inv(P(:,:,cnts(tokeep(jj))));
        xhat = x(:,cnts(tokeep(jj)));
       hatMatrix = hatMatrix+tau(jj).*[PiHat, -PiHat*xhat; -xhat'*PiHat, xhat'*PiHat*xhat+c(tokeep(jj))];
    end

    LMI = tau>=0;
    LMI = [LMI,barMatrix-hatMatrix>=0];
    LMI = [LMI,sum(tau)==1];
    d1 = optimize(LMI,[],opts);
%     value(tau)
    
    if(d1.problem~=0)
%     if(true)
        tokeep = [tokeep,kk];
        Pset_tmp = [Pset_tmp,[Pset_prev(:,hst(kk));nxt(kk)]];
    end
    
end

xout = x(:,cnts(tokeep));
Pout = P(:,:,cnts(tokeep));
cout = c(tokeep);
Pset_out = Pset_tmp;
ind_min = nxt(1);
end