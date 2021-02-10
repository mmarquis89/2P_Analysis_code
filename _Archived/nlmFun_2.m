function [yhat, test] = nlmFun_2(beta, x)

% % Fit tau
% b1 = beta(1); % f
% b2 = beta(2); % tau (in volumes, ~160 msec/vol)
% b3 = beta(3); % adaptation term coefficient
% b4 = beta(4); % intercept term
% 
% Fixed tau
b1 = beta(1); % f
b2 = beta(2); % adaptation term coefficient
b3 = beta(3); % intercept

x1 = x(:,1);
x2 = x(:,2);

% Return result that will have high error if parameters are outside acceptable ranges
if b1 > 1 || b1 < 0
    yhat = ones(size(x1)) * 1000;
    return;
end

% % Fit tau
% if b2 < 1 || b2 > numel(x1)*1
%     yhat = ones(size(x1)) * 1000;
%     return;
% end

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
%     fcn = @(x) (1 - x) ./ b2; % Fit tau    
    fcn = @(x) (1 - x) ./ 3000; % Fixed tau
    currA = currA + fcn(currA);
    yhat(end + 1) = x1(iVol) .* currA;
end
% yhat = b3*yhat + b4; % Fit tau
yhat = b2*yhat + b3;   % Fixed tau
yhat = yhat';

end