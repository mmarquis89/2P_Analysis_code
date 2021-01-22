function f = makefig(figNum)
% Simply creates a figure (or clears it if it already exists), sets "hold on", changes the 
% background color to white, and returns the handle. Will call "figure(1)" with no argument, or you 
% can provide a different figure number.
if nargin < 1
    figNum = 1;
end

f = figure(figNum);
f.Color = [1 1 1];
clf;
hold on;


end