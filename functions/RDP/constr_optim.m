function [c,ceq] = constr_optim(x,xhat,eta)
    c = norm(x-xhat,2)-eta;
    ceq = [];
end