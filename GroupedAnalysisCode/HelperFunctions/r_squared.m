function [r2,r2adj] = r_squared(ydata,yestimation,nparam)
% RSQUARED calculates the coefficient of determination (r2) from the original
%   data (ydata) and fited data (yestimation). It also calculates the adjusted
%   coefficient (r2adj) considering the number of parameters of the model
%   (nparam).
%
% Syntax:  
%    r2 = rsquared(ydata,yestimation)
%    [r2,r2adj]=rsquared(ydata,yestimation,nparam)
%
% Example:  
%    xdata = [1    5  14  23  25  48   49   59   73   77   99  ];
%    ydata = [-100 70 100 450 550 2200 2300 3400 5300 5906 9600];
%    plot(xdata,ydata,'ok'), hold on
%    param_1 = polyfit(xdata,ydata,1);
%    yestimation_1 = polyval(param_1,xdata);
%    [r2_1,r2adj_1]=rsquared(ydata,yestimation_1,length(param_1))
%    plot(xdata,yestimation_1,'-r')
%    param_2 = polyfit(xdata,ydata,2);
%    yestimation_2 = polyval(param_2,xdata);
%    plot(xdata,yestimation_2,'-b')
%    [r2_2,r2adj_2]=rsquared(ydata,yestimation_2,length(param_2))
%    legend({'data',['r2=' num2str(r2_1) ', r2adj=' ...
%        num2str(r2adj_1)],['r2=' num2str(r2_2) ', r2adj=' num2str(r2adj_2)]}, ...
%        'Location','best')
%
% Equations
%    SSres=sum( (ydata-yestimation).^2 );                                 % residual sum of squares
%    SStot=sum( (ydata-mean(ydata)).^2 );                                 % total sum of squares
%    r2=1-SSres/SStot;                                                    % standard rsquared
%   r2adj = 1 - SSres/SStot * (length(ydata)-1)/(length(ydata)-nparam);   % adjust for the number of parameters
%
% Check https://en.wikipedia.org/wiki/Coefficient_of_determination
%
% Comments? rpavao@gmail.com
if nargin<2
    help rsquared
else
    SSres=sum( (ydata-yestimation).^2 );
    SStot=sum( (ydata-mean(ydata)).^2 );
    %SStot = (length(ydata)-1) * var(ydata);
    
    r2=1-SSres/SStot;
    
    if nargout>1
        if nargin<3
            warning('using nparam=2 for caclulating r2adj (ex. 1st-degree polynomial fit)')
            nparam=2;
        end
        r2adj = 1 - SSres/SStot * (length(ydata)-1)/(length(ydata)-nparam);        
    end
end

% % 
% % calcuateR2 Cacluate R-squared
% % R2 = calcuateR2(z,z_est) takes two inputs - The observed data x and its
% % estimation z_est (may be from a regression or other model), and then
% % compute the R-squared value a measure of goodness of fit. R2 = 0
% % corresponds to the worst fit, whereas R2 = 1 corresponds to the best fit.
% % 
% % Copyright @ Md Shoaibur Rahman (shaoibur@bcm.edu)
% r = z-zEst;
% normr = norm(r);
% SSE = normr.^2;
% SST = norm(z-mean(z))^2;
% R2 = 1 - SSE/SST;
end