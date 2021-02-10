function yhat = nlmFun_odorOnly(beta, x)

% Read in starting values for fit parameters
b_odor = beta(1);        % odor term coefficient
intercept = beta(2);     % intercept term

% Get predictor variables
odorResp = x(:, 1);    % Estimated odor responses

% Generate final output prediction
yhat = b_odor*odorResp' + intercept;
yhat = yhat';

end