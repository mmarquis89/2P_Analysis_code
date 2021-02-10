function yhat = nlmFun(beta, x)

% fl ~ b1*moveSpeed + b2*odorResp + (b3*odorResp*(odorHist_120^b4))
% y ~ b1*x1 + b2*x2 + (b3*x3*(x4^b4))

b1 = beta(1);
b2 = beta(2);
b3 = beta(3);
b4 = beta(4);
b5 = beta(5);
b6 = beta(6);

x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);

if b5 < 0
    yhat = ones(size(x, 1), 1) * 100;
else
%     yhat =  b1.*x1 + b2.*x2 + b3.*x3 + (b4.*x2.*(x3.^b5)) + b6;
%     yhat =  b1.*x1 + b2.*x2 + b3.*x3 + (b4.*x2.*(1-exp(b5.*x3))) + b6; 
    fcn = @(x) x ./ (1 + x.^b5).^(1./b5);
    yhat =  b1.*x1 + b2.*x2 + b3.*x3 + (b4.*x2.*fcn(x3)) + b6;
end

end