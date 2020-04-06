
parentDir = 'C:\Users\Wilson Lab\Documents\MATLAB';
fileName = fullfile(parentDir, 'timeSeries.csv');

% Download newest version of file
test = websave(fileName, ['https://raw.githubusercontent.com/', ...
        'CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_', ...
        'covid19_confirmed_global.csv']);

% Read and format table
tb = readtable(fileName, 'ReadVariableNames', 1);
tb = tb(:, [1 2 5:size(tb, 2)]);
tb.Properties.VariableNames{1} = 'Province';
tb.Properties.VariableNames{2} = 'Country';
monthStr = {'Jan', 'Feb', 'Mar', 'Apr', 'May'};
for iCol = 3:size(tb, 2)
    currName = tb.Properties.VariableNames{iCol};
    currMonthName = regexp(currName, '(?<=x)..?(?=_)', 'match');
    currMonthName = monthStr{str2double(currMonthName{:})};
    moName = regexprep(currName, '(?<=x)..?(?=_)', currMonthName);
    newName = regexp(moName(2:end), '..._..?(?=_)', 'match');
    tb.Properties.VariableNames(iCol) = newName;
end


%% Extract data for target countries 

countryList = {'US', 'Italy', 'Spain', 'Korea, South', 'Iran', 'France', 'Germany', 'United Kingdom'};
% countryList = {'Brazil', 'Switzerland', 'United Kingdom', 'Russia'};
startDate = 'Mar_2';

for iCountry = 1:numel(countryList)
    currName = countryList{iCountry};
    currTb = tb(strcmp(tb.Country, currName), :);
    currTbArr = [currName, num2cell(sum(currTb{:, 3:end}, 1))];
       
    if iCountry == 1
        outputTb = cell2table(currTbArr, 'VariableNames', currTb.Properties.VariableNames(2:end));
    else
        outputTb = [outputTb; cell2table(currTbArr, 'VariableNames', ...
                currTb.Properties.VariableNames(2:end))];
    end
    
end

% Plot data
f = figure(1); clf;
f.Color = [1 1 1];
plotData = outputTb{:, 2:end}';
plot(plotData, '-o', 'linewidth', 1);
plotDates = outputTb.Properties.VariableNames(2:end);
if isempty(startDate)
    xlim([0, size(plotData, 1) + 2]);
else
    xlim([find(strcmp(plotDates, startDate)) - 1, size(plotData, 1) + 2]);
end
ylim([-0.03, 1.05] * max(plotData(:)));
legend(outputTb.Country, 'location', 'nw')
ax = gca;
ax.XTick = 1:2:size(plotData, 1);
ax.XTickLabel = regexprep(plotDates(1:2:end), '_', '-');
ax.XTickLabelRotation = -60;
ax.YTickLabel = ax.YTick;

test = diff(plotData, 1);% ./ plotData(1:end-1, :);
figure(2);clf; plot(smoothdata(test, 1, 'gaussian', 3), '-o')
xlim([find(strcmp(plotDates, startDate)) - 1, size(plotData, 1) + 2]);

