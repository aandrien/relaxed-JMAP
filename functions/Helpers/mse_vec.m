function y = mse_vec(x,u)
% Mean squared error (MSE) for vectors
y = (norm(x-u,2)^2)/length(x);
end