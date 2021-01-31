function [yhat, yhat_odor, yhat_moveSpeed] = nlmFun_double_tau_fit(beta, x)

% Read in starting values for fit parameters
f_odor = beta(1);        % f_odor
f_moveSpeed = beta(2);   % f_moveSpeed

tau_odor = beta(3);      % tau_odor (in volumes, ~160 msec/vol)
tau_moveSpeed = beta(4); % tau_moveSpeed

b_odor = beta(5);        % odor adaptation term coefficient
b_moveSpeed = beta(6);   % odor adaptation term coefficient

b_odorHistory = beta(7);    % odor history term coefficient

intercept = beta(8);     % intercept term

% Get predictor variables
odorResp = x(:, 1);    % Estimated odor responses
moveSpeed = x(:, 2);   % moveSpeed
odorOnsets = x(:, 3);  % odor onset vector
odorHistory = x(:, 4); % odor history

% Return result that will have high error if any parameters are outside acceptable ranges
if f_odor > 1 || f_odor < 0 || f_moveSpeed > 1 || f_moveSpeed < 0
    yhat = ones(size(odorResp)) * 1000;
    return;
end
if tau_odor < 1 || tau_odor > numel(odorResp)*1 || tau_moveSpeed < 1 || tau_moveSpeed > numel(odorResp)*1
    yhat = ones(size(odorResp)) * 1000;
    return;
end
if b_odorHistory > 0
    yhat = ones(size(odorResp)) * 1000;
    return;
end

% Calculate odor response w/adaptation throughout experiment
firstOdor = 0;
currA_odor = 1;
currA_moveSpeed = 1;
yhat_odor = [];
yhat_moveSpeed = [];
for iVol = 1:numel(odorResp)
    
    % Apply adaptation terms
    if odorOnsets(iVol)
        if firstOdor
            % Only start decreasing the response amplitude after the first odor stim
            currA_odor = currA_odor .* f_odor;
            currA_moveSpeed = currA_moveSpeed .* f_moveSpeed;
        else
            firstOdor = 1;
        end
    end
    
    % Apply recovery terms
    fcn_odor = @(x) (1 - x) ./ tau_odor;
    currA_odor = currA_odor + fcn_odor(currA_odor);
    fcn_moveSpeed = @(x) (1 - x) ./ tau_moveSpeed;
    currA_moveSpeed = currA_moveSpeed + fcn_moveSpeed(currA_moveSpeed);
    
    % Calculate adapted responses
    yhat_odor(end + 1) = odorResp(iVol) .* currA_odor;
    yhat_moveSpeed(end + 1) = moveSpeed(iVol) .* currA_moveSpeed;
end

% Generate final output prediction
yhat = b_odor*yhat_odor + b_moveSpeed*yhat_moveSpeed + b_odorHistory*odorHistory' + intercept;
yhat = yhat';

end