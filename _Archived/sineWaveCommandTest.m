


s = daq.createSession('ni');
s.Rate = 100000;
s.addAnalogOutputChannel('Dev1', 0, 'Voltage');
s.addAnalogInputChannel('Dev1', 0, 'Voltage');
duration = 5;
t = 0:(1/s.Rate):duration;
f = 80;
x = sin(2*pi*t*f);

outputData = x * 10;

s.queueOutputData(outputData')
test = s.startForeground();


