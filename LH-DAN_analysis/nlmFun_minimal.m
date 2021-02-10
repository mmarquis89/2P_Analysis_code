function yhat = nlmFun_minimal(beta, x)

% Read in starting values for fit parameters
b_odor = beta(1);        % odor term coefficient
b_moveSpeed = beta(2);   % odor term coefficient
intercept = beta(3);     % intercept term

% Get predictor variables
odorResp = x(:, 1);    % Estimated odor responses
moveSpeed = x(:, 2);   % moveSpeed

% Generate final output prediction
yhat = b_odor*odorResp' + b_moveSpeed*moveSpeed' + intercept;
yhat = yhat';

end