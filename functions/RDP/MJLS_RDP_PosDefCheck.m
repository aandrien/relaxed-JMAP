function [xout,Pout,cout,ind_min,Pset_out] = MJLS_RDP_PosDefCheck(sys,xprev,Pprev,cprev,yk1,uk,uk1,Pset_prev,e)

n = sys.n;
s = sys.s;

%% Calculate KF for each entry in set C
cnt = 1;
hst = zeros(size(Pset_prev,2)*s,1);
nxt = zeros(size(Pset_prev,2)*s,1);
x = zeros(n,size(Pset_prev,2)*s);
P = zeros(n,n,size(Pset_prev,2)*s);
c = zeros(size(Pset_prev,2)*s,1);

for ii = 1:size(Pset_prev,2) %Loop over entries in P_{k-1}
    for jj = 1:s %For each entry in P_{k-1} there are s new options (set C_K)
        [x(:,cnt),P(:,:,cnt),c(cnt)] = MJLS_KF2(sys,xprev(:,ii),Pprev(:,:,ii),Pset_prev(end,ii),jj,yk1,uk,uk1,cprev(ii),0);
        hst(cnt) = ii;
        nxt(cnt) = jj;
        cnt = cnt+1;
    end
end
cnt = cnt-1;

%% RDP procedure using minmax
%[1] initialize set Pset as empty
tokeep = zeros(1,cnt); %We don't know how big Pset will be, initialize to max possible for speed
hatMatrices = zeros(n+1,n+1,cnt);

%[2] Add the element in Cset with the smallest c to Pset
c_hst_nxt = sortrows([c,hst,nxt,[1:cnt]'],1); %sort
c = c_hst_nxt(:,1);
hst = c_hst_nxt(:,2);
nxt = c_hst_nxt(:,3);
cnts = c_hst_nxt(:,4);

tokeep(1) = 1;
nP = 1;

PiHat = inv(P(:,:,cnts(tokeep(nP))));
xhat = x(:,cnts(tokeep(nP)));
% Next two lines are necessary to make sure the hatMatrix is
% symmetric, since else yalmip will consider it to be full.
Qtmp = [PiHat./2, -PiHat*xhat./2; -xhat'*PiHat./2, xhat'*PiHat*xhat./2+c(tokeep(nP))];
% Qtmp = [PiHat./2, -PiHat*xhat./2; -xhat'*PiHat./2, xhat'*PiHat*xhat./2+abs(c(tokeep(nP)))];
Qtmp =(triu(Qtmp)+triu(Qtmp)') - eye(size(triu(Qtmp),1)).*diag(triu(Qtmp));
% Add this entry of Pset to hatMatrices
hatMatrices(:,:,nP) = Qtmp;

%[3] Take one element i in Ck and check if it satisfies relaxed cost
for kk = 2:cnt
    %Calculate relaxed cost
    [xbar,Pbar,cbar] = MJLS_KF2(sys,xprev(:,hst(kk)),Pprev(:,:,hst(kk)),Pset_prev(end,hst(kk)),nxt(kk),yk1,uk,uk1,cprev(hst(kk)),e);
    
    %Check condition by comparing eigenvalues of matrices
    PiBar = inv(Pbar);
    barMatrix = [PiBar./2, -PiBar*xbar./2; -xbar'*PiBar./2, xbar'*PiBar*xbar./2+cbar];
%     barMatrix = [PiBar./2, -PiBar*xbar./2; -xbar'*PiBar./2, xbar'*PiBar*xbar./2+abs(cbar)];
    barMatrix =(triu(barMatrix)+triu(barMatrix)') - eye(n+1).*diag(triu(barMatrix)); % Make sure barMatrix is symmetric!
    if(~issymmetric(barMatrix))
        disp('Error: barMatrix not symmetric')
    end
    
    valid = false;
    for ii = 1:nP
        if(all(eig(barMatrix-hatMatrices(:,:,ii))>=0))
            valid = true;
            break
        end
    end
    
    
    if(~valid) %Condition is not met, so add to set P
        nP = nP+1;
        tokeep(nP) = kk;
        
        PiHat = inv(P(:,:,cnts(tokeep(nP))));
        xhat = x(:,cnts(tokeep(nP)));
        % Next two lines are necessary to make sure the hatMatrix is
        % symmetric, since else yalmip will consider it to be full.
%         Qtmp = [PiHat./2, -PiHat*xhat./2; -xhat'*PiHat./2, xhat'*PiHat*xhat./2+c(tokeep(nP))];
        Qtmp = [PiHat./2, -PiHat*xhat./2; -xhat'*PiHat./2, xhat'*PiHat*xhat./2+c(tokeep(nP))];
        Qtmp =(triu(Qtmp)+triu(Qtmp)') - eye(size(triu(Qtmp),1)).*diag(triu(Qtmp));
        if(~issymmetric(Qtmp))
            disp('Error: hatMatrix not symmetric')
        end
        % Add this entry of Pset to hatMatrices
        hatMatrices(:,:,nP) = Qtmp;
    end
    
end
tokeep = tokeep(1:nP); %remove zero entries from tokeep

xout = x(:,cnts(tokeep));
Pout = P(:,:,cnts(tokeep));
cout = c(tokeep);
Pset_out = [Pset_prev(:,hst(tokeep));nxt(tokeep)'];
ind_min = nxt(1);
% if(size(tokeep,2)>=4)
%     plotquadcost(x,P,c)
% end
end