function yhat = nlmFun_no_adaptation(beta, x)

% Read in starting values for fit parameters
b_odor = beta(1);        % odor term coefficient
b_moveSpeed = beta(2);   % odor term coefficient
b_odorHistory = beta(3); % odor history term coefficient
intercept = beta(4);     % intercept term

% Get predictor variables
odorResp = x(:, 1);    % Estimated odor responses
moveSpeed = x(:, 2);   % moveSpeed
odorHistory = x(:, 3); % odor history

% Return result that will have high error if any parameters are outside acceptable ranges
if b_odorHistory > 0
    yhat = ones(size(odorResp)) * 1000;
    return;
end

% Generate final output prediction
yhat = b_odor*odorResp' + b_moveSpeed*moveSpeed' + b_odorHistory*odorHistory' + intercept;
yhat = yhat';

end