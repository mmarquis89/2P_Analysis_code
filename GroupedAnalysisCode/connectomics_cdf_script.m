dt = DataTable(tbl);
filterDefs = struct();
typeList = unique(dt.sourceData.neuronType);


currTypeList = typeList([1:7, 11:13, 8:10]);
partnerDirection = 'upstream';
yVarType = 'count'; % 'percentage'; %
yScaleType =  'log'; %'linear'; % 
xNorm = 1;
xMinPct = 70;
groupPartnersByType = 0;
tracedPartnersOnly = 1;
legendPos = 'nw';
saveFig = 0;

cm = jet(numel(currTypeList)) * 0.9;

% Create figure
f = figure(1);clf;
f.Color = [1 1 1];
ax = axes(); hold on;

filterDefs.partnerDirection = partnerDirection;
if strcmp(filterDefs.partnerDirection, 'downstream')
    directionStr = 'output';
else
    directionStr = 'input';
end

partnerTypeCounts = [];
maxSynapseCounts = [];
totalSynapseCounts = [];

overallMinVal = inf;
overallMaxVal = 0;
for iType = 1:numel(currTypeList)
    
    filterDefs.neuronType = currTypeList{iType};
    
    % Get data subset
    dt = dt.initialize_filters(filterDefs);
    sub = dt.subset;
    
    % Drop rows for partners that are incompletely traced or that do not have an assigned cell type
    if tracedPartnersOnly
        sub = sub(strcmpi(sub.status, 'traced') & ~cellfun(@isempty, sub.partnerType) & ...
            ~strcmp(sub.partnerType, 'NA'), :);
    end
        
    % Calculate connection weight of with each partner neuron as a percentage of total DAN synapses
    totalSynapses = sum(sub.synapseCount);
    maxTotalSynapses(iType) = max(totalSynapses);
    connectionWeightPct = 100 * (sub.synapseCount ./ totalSynapses);
    sub.connectionWeightPct = connectionWeightPct;
    
    % Create another table that groups partners by cell type
    groupedSub = groupsummary(sub, {'neuronID', 'neuronName', 'partnerDirection', 'partnerType'}, ...
            'sum', {'synapseCount'});
    groupedSub = groupedSub(~strcmp(groupedSub.partnerType, 'NA'), :);
    groupedSub.connectionWeightPct = 100 * (groupedSub.sum_synapseCount ./ totalSynapses);
    
    totalSynapseCounts(iType) = sum(groupedSub.sum_synapseCount);
    partnerTypeCounts(iType) = numel(groupedSub.sum_synapseCount);
    maxSynapseCounts(iType) = max(groupedSub.sum_synapseCount);
    
    % Get correct plot data for selected options
    if groupPartnersByType
        pctData = groupedSub.connectionWeightPct;
        countData = groupedSub.sum_synapseCount;
    else
        pctData = sub.connectionWeightPct;
        countData = sub.synapseCount;
    end
    if strcmpi(yVarType, 'count')
        yy = sort(countData);
    else
        yy = sort(pctData);
    end
    if xNorm
        xx = (1:numel(yy)) ./ numel(yy);
    else
        xx = 1:numel(yy);
    end
    
    % Keep track of overall min and max vals for setting final axes limits
    overallMinVal = min(overallMinVal, min(yy));
    overallMaxVal = max(overallMaxVal, max(yy));
    
    plot(ax, xx, yy, '-', 'linewidth', 1.5, 'color', cm(iType, :));
    ax.Children(1).Marker = '.';
    ax.Children(1).MarkerSize = 12;
    ax.FontSize = 14; 
     
    if iType == 1
        lineObjects = ax.Children;
    else
        lineObjects(end + 1) = ax.Children(1);
    end
    
end

% Change Y-axis to a log scale if necessary
if strcmp(yScaleType, 'log')
    ax = logify(ax, [log10(exp_floor(overallMinVal)), log10(exp_ceil(overallMaxVal))]);
end

% Set Y limits
if strcmp(yScaleType, 'log')
    ax.YLim = [overallMinVal, overallMaxVal] .* [0.9, 1.5];
else
    ax.YLim = [overallMinVal, overallMaxVal] .* [0.9, 1.1];
end

% Set X limits if necessary
if xNorm
    xlim(ax, [xMinPct/100, 1.02]);
%     ax.XTickLabels = '';
elseif xMinPct > 0
    xL = xlim();
    pctThresh = round(xL(2) * (xMinPct / 100));
    xlim([pctThresh, xL(2)])
end

% Add legend and labels
if strcmp(yVarType, 'count')
    ylabel('Synapse count');
else
    ylabel(['% of total ', directionStr, ' synapses']);
end
if groupPartnersByType
    xlabel([capitalize(partnerDirection), ' partner cell types'])
else
    xlabel([capitalize(partnerDirection), ' partner cells'])
end
if xNorm
    ax.XLabel.String = [ax.XLabel.String, ' (percentile)'];
end
legend(lineObjects, regexprep(currTypeList, '_', '\\_'), 'location', legendPos)
title([capitalize(partnerDirection), ' synaptic partner connection strengths'])

if saveFig
    saveDir = fullfile(parentDir, 'figs');
    
    if groupPartnersByType
        partnerTypeStr = 'cellTypes';
    else
        partnerTypeStr = 'cells';
    end
    if xNorm
        xNormStr = 'normalized';
    else
        xNormStr = 'raw';
    end
    
    fileName = [partnerDirection, '_partner_', partnerTypeStr, '_synapse', capitalize(yVarType), ...
        's_', yScaleType, 'Y_', xNormStr, 'X_xMin_', num2str(xMinPct)];
    
    save_figure(f, saveDir, fileName);
end


% Create summary tables of the top partner cells and cell types
disp(' '); disp(' ')

nRows = 20;

cellTypeTbl = tail(sortrows(groupedSub, 'connectionWeightPct'), nRows);
cellTypeTbl.Properties.VariableNames{strcmp(cellTypeTbl.Properties.VariableNames, ...
        'sum_synapseCount')} = 'synapseCount';
cellTypeTbl.Properties.VariableNames{strcmp(cellTypeTbl.Properties.VariableNames, ...
        'GroupCount')} = 'nCells';
    
singleCellTbl = tail(sortrows(sub(:, {'neuronID', 'neuronName', 'partnerName', 'partnerType', ...
        'partnerDirection', 'synapseCount', 'connectionWeightPct'}), 'connectionWeightPct'), nRows);
    
disp(cellTypeTbl)
disp(singleCellTbl);


%% Create scatter plots of summary variables
tb = table(currTypeList, totalSynapseCounts', partnerTypeCounts', maxSynapseCounts', ...
        'variablenames', {'cellType', 'totalSynapses', 'totalPartnerTypes', 'maxSynapseCount'});

xx = tb.totalSynapses;
yy = tb.maxSynapseCount;
xLab = ['Total ', directionStr, ' synapses'];
yLab = ['Max synapses with a single partner cell type'];
saveFig = 1;
fileName = ['total_', directionStr, '_synapses_vs_max_synapses_with_single_cell_type'];

dx = ones(size(xx)) * 100;
% dx(7) = dx(7) - 850;
dy = ones(size(xx)) * 25;
dy(9) = dy(9) - 50;

f = figure(4);clf;
f.Color = [1 1 1];
hold on; 
plot(xx, yy, '.', 'color', 'b', 'markersize', 25);  
ax = gca();
ax.FontSize = 14;
xlabel(xLab)
ylabel(yLab)
c = regexprep(tb.cellType, '_', '\\_');
text(xx + dx, yy + dy, c, 'FontSize', 12);

if saveFig
    save_figure(f, saveDir, fileName);
end

%% Local functions 

% Change Y axis to a log scale and set tick labels
function ax = logify(ax, yExpLims)

YTickVals = [];
baseYTickVals = [10^yExpLims(1):10^yExpLims(1):(10^yExpLims(1) * 5)];
for i = 1:diff(yExpLims)
    YTickVals = [YTickVals, baseYTickVals * 10^(i-1)];
end
YTickVals(end + 1) = 10^yExpLims(2);
ax.YScale = 'log';
ax.YTick = YTickVals;

end

% Round up to the next power of 10
function x = exp_ceil(n, startingExp)
    if nargin < 2
        startingExp = -5;
    end
    if n > 10^startingExp
        x = exp_ceil(n, startingExp + 1);        
    else
        x = 10^startingExp;
    end
end

% Round down to the next power of 10
function x = exp_floor(n, startingExp)
    if nargin < 2
        startingExp = 10;
    end
    if n < 10^startingExp
        x = exp_floor(n, startingExp - 1);        
    else
        x = 10^startingExp;
    end
end
