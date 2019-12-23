function [xout,Pout,cout,ind_min,Pset_out] = MJLS_RDP_minmax(xprev,Pprev,cprev,yk1,uk,uk1,Pset_prev,e,eta)

global s

n = size(xprev,1);
%% Calculate KF for each entry in set C
cnt = 1;
hst = zeros(size(Pset_prev,2)*s,1);
nxt = zeros(size(Pset_prev,2)*s,1);
x = zeros(n,size(Pset_prev,2)*s);
P = zeros(n,n,size(Pset_prev,2)*s);
c = zeros(size(Pset_prev,2)*s,1);

for ii = 1:size(Pset_prev,2) %Loop over entries in P_{k-1}
    for jj = 1:s %For each entry in P_{k-1} there are s new options (set C_K)
        [x(:,cnt),P(:,:,cnt),c(cnt)] = MJLS_KF(xprev(:,ii),Pprev(:,:,ii),Pset_prev(end,ii),jj,yk1,uk,uk1,cprev(ii),0);
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
xhat = x(:,cnts(tokeep));
Pset_tmp = [Pset_prev(:,hst(1));nxt(1)];

%[3] Take one element i in Ck and check if it satisfies relaxed cost
options = optimoptions('fminimax','Display','off',...
    'StepTolerance',1e-12,...
    'FunctionTolerance',1e-12,...
    'OutputFcn',@outfun);
% 'MaxFunctionEvaluations',50000,...
%     'MaxIterations',5000,...
for kk = 2:cnt-1
    %Calculate relaxed cost
    [xbar,Pbar,cbar] = MJLS_KF(xprev(:,hst(kk)),Pprev(:,:,hst(kk)),Pset_prev(end,hst(kk)),nxt(kk),yk1,uk,uk1,cprev(hst(kk)),e);
    
    %Check condition
    [~,~,theta,exflag] = fminimax(@(xv) func_optim(xv,xbar,Pbar,cbar,x(:,cnts(tokeep)),P(:,:,cnts(tokeep)),c(tokeep)),xhat,[],[],[],[],[],[],@(xv) constr_optim(xv,xhat,eta),options);
    %     exflag  %      1  fminimax converged to a solution.
    %       4  Computed search direction too small.
    %       5  Predicted change in max objective function too small.
    %       0  Too many function evaluations or iterations.
    %      -1  Stopped by output/plot function.
    %      -2  No feasible point found.
    if(theta <= 0 || exflag == 0 || exflag == -2 || exflag == -1)
        %     if(true)
        tokeep = [tokeep,kk];
        Pset_tmp = [Pset_tmp,[Pset_prev(:,hst(kk));nxt(kk)]];
    end
%     if(exflag==1 || exflag == -1)
%         theta
%         exflag
%     end
end

xout = x(:,cnts(tokeep));
Pout = P(:,:,cnts(tokeep));
cout = c(tokeep);
Pset_out = Pset_tmp;
ind_min = nxt(1);
end