function F = func_optim(x,xbar,Pbar,cbar,xhat,Phat,chat)
    sizeP = size(xhat,2);
    F = zeros(sizeP,1);
    for ii = 1:sizeP
       F(ii) = 0.5*mnorm(x-xbar,inv(Pbar))^2+cbar-0.5*mnorm(x-xhat(:,ii),inv(Phat(:,:,ii)))^2-chat(ii);
    end
end