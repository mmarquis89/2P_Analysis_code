function yhat = nlmFun_speedOnly(beta, x)

% Read in starting values for fit parameters
b_moveSpeed = beta(1);   % speed term coefficient
intercept = beta(2);     % intercept term

% Get predictor variables
moveSpeed = x(:, 1);    % Estimated odor responses

% Generate final output prediction
yhat = b_moveSpeed*moveSpeed' + intercept;
yhat = yhat';

end