function [xout,Pout,cout,ind_min,Pset_out] = MJLS_RDP_LMI(xprev,Pprev,cprev,yk1,uk,uk1,Pset_prev,e)

global s

n = size(xprev,1);

%% SOLVER SETTINGS
% opts = sdpsettings('verbose',0,'solver','sdpt3');
% opts =sdpsettings('verbose',0,'solver', 'sedumi', 'sedumi.eps', 1e-8, ...
%                 'sedumi.cg.qprec', 1, 'sedumi.cg.maxiter', 49, ...
%                 'sedumi.stepdif', 2);
opts = sdpsettings('verbose',0,'solver','mosek');

%% Calculate KF for each entry in set C
cnt = 1;
hst = zeros(size(Pset_prev,2)*s,1);
nxt = zeros(size(Pset_prev,2)*s,1);
x = zeros(n,size(Pset_prev,2)*s);
P = zeros(n,n,size(Pset_prev,2)*s);
c = zeros(size(Pset_prev,2)*s,1);

for ii = 1:size(Pset_prev,2) %Loop over entries in P_{k-1}
    for jj = 1:s %For each entry in P_{k-1} there are s new options (set C_K)
        [x(:,cnt),P(:,:,cnt),c(cnt)] = MJLS_KF2(xprev(:,ii),Pprev(:,:,ii),Pset_prev(end,ii),jj,yk1,uk,uk1,cprev(ii),0);
        hst(cnt) = ii;
        nxt(cnt) = jj;
        cnt = cnt+1;
    end
end
cnt = cnt-1;

%% RDP procedure using minmax
%[1] initialize set Pset as empty
tokeep = zeros(1,cnt); %We don't know how big Pset will be, initialize to max possible for speed

%[2] Add the element in Cset with the smallest c to Pset
c_hst_nxt = sortrows([c,hst,nxt,[1:cnt]'],1); %sort
c = c_hst_nxt(:,1);
hst = c_hst_nxt(:,2);
nxt = c_hst_nxt(:,3);
cnts = c_hst_nxt(:,4);

tokeep(1) = 1;
nP = 1;

%[3] Take one element i in Ck and check if it satisfies relaxed cost
for kk = 2:cnt
    %Calculate relaxed cost
    [xbar,Pbar,cbar] = MJLS_KF2(xprev(:,hst(kk)),Pprev(:,:,hst(kk)),Pset_prev(end,hst(kk)),nxt(kk),yk1,uk,uk1,cprev(hst(kk)),e);
    
    %Check condition using LMIs
    tau = sdpvar(1,nP,'full');
    PiBar = inv(Pbar);
    barMatrix = [PiBar./2, -PiBar*xbar./2; -xbar'*PiBar./2, xbar'*PiBar*xbar./2+cbar];
    barMatrix =(triu(barMatrix)+triu(barMatrix)') - eye(n+1).*diag(triu(barMatrix)); % Make sure barMatrix is symmetric!
    if(~issymmetric(barMatrix))
        disp('Error: barMatrix not symmetric')
    end
    
    hatMatrix = zeros(n+1,n+1);
    for jj = 1:nP
        PiHat = inv(P(:,:,cnts(tokeep(jj))));
        xhat = x(:,cnts(tokeep(jj)));
        % Next two lines are necessary to make sure the hatMatrix is
        % symmetric, since else yalmip will consider it to be full.
        Qtmp = [PiHat./2, -PiHat*xhat./2; -xhat'*PiHat./2, xhat'*PiHat*xhat./2+c(tokeep(jj))];
        Qtmp =(triu(Qtmp)+triu(Qtmp)') - eye(size(triu(Qtmp),1)).*diag(triu(Qtmp));
        hatMatrix = hatMatrix+tau(jj)*Qtmp;
        
        if(~issymmetric(hatMatrix))
            disp('Error: hatMatrix not symmetric')
        end
        
    end
    
    LMI = tau>=0;
    LMI = [LMI,barMatrix-hatMatrix>=0];
    LMI = [LMI,sum(tau)==1];
    d1 = optimize(LMI,[],opts);
    
    if(d1.problem~=0) %Condition is not met, so add to set P
        nP = nP+1;
        tokeep(nP) = kk;
    end
    if(d1.problem ~= 0 && d1.problem~=1)
        d1
    end
end
tokeep = tokeep(1:nP); %remove zero entries from tokeep

xout = x(:,cnts(tokeep));
Pout = P(:,:,cnts(tokeep));
cout = c(tokeep);
Pset_out = [Pset_prev(:,hst(tokeep));nxt(tokeep)'];
ind_min = nxt(1);
end