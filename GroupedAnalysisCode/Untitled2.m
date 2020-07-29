flt = [1 1 1 1 1];
input = [1 1 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 1 0 1 1 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0];
% flt = ones(size(input));

figure(10);clf; hold on;

sd = 1
mu = 0
x = linspace(-4*sd,4*sd,10);
y = 1/(2*pi*sd)*exp(-(x-mu).^2/(2*sd^2));
flt = y;
plot(flt)


plot(input, '-o'); 
ylim([0 5])

output = conv(flt, input);
plot(output, '-o')



% input = testOdorData(1:500);
input = [ones(1, 14), zeros(1, 28)];
output = testFl;

[q, r] = deconv(output, input(1:end-1) + 1);

sd = 1
mu = 0
x = linspace(-4*sd,4*sd,10);
y = 1/(2*pi*sd)*exp(-(x-mu).^2/(2*sd^2));
% plot(y);

testGaussian = y;










