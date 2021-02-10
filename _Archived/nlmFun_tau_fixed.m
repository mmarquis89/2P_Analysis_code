function yhat = nlmFun_tau_fixed(beta, x)

% Set tau
FIXED_TAU = 1000;

b1 = beta(1); % f
b2 = beta(2); % adaptation term coefficient
b3 = beta(3); % intercept

x1 = x(:,1); % Estimated odor responses
x2 = x(:,2); % Odor onset vector

% Return result that will have high error if parameters are outside acceptable ranges
if b1 > 1 || b1 < 0
    yhat = ones(size(x1)) * 1000;
    return;
end

% Calculate odor response w/adaptation throughout experiment
firstOdor = 0;
currA = 1;
yhat = [];
for iVol = 1:numel(x1)
    if x2(iVol)
        if firstOdor
            % Only start decreasing the response amplitude after the first odor stim
            currA = currA .* b1;
        else
            firstOdor = 1;
        end
    end
    fcn = @(x) (1 - x) ./ FIXED_TAU;
    currA = currA + fcn(currA);
    yhat(end + 1) = x1(iVol) .* currA;
end
yhat = b2*yhat + b3; 
yhat = yhat';

end