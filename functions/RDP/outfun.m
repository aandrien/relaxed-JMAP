function stop = outfun(x,optimValues,state)
% Stop searching as soon as function value below zero is found that
% satisfies constraints.
if(max(optimValues.fval)<=0 && strcmp(state, 'iter') && optimValues.constrviolation<=0)
    stop = true;
else
    stop = false;
end

end