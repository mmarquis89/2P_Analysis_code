
flowRate = 20; % mL/min
tubeDia = g20; % in



% Convert flow to cubic in/sec
flowRate = flowRate / 983.224;

% Multiply by area at tube mouth to get flow speed in in/sec
tubeArea = pi * ((tubeDia / 2)^2);

% Calculate flow speed in cm/sec
flowSpeed = (flowRate / tubeArea) * 2.54;

disp(' ')
disp(['Flow rate: ', num2str(flowRate, 3)])
disp(['Tube diameter: ', num2str(tubeDia, 3)])
disp(['Tube area: ', num2str(tubeArea, 3)])
disp(['Flow speed: ', num2str(flowSpeed, 3)])