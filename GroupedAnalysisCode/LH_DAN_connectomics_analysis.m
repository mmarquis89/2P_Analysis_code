
% CONCATENATE AND RE-SAVE CONNECTIVTY DATA FOR ALL LH-DANS

% Load all connectivity data
parentDir = 'C:\Users\Wilson Lab\Google Drive\Lab Work\LH-DAN_connectomics_analysis';
neuron_names = {'TypeB', 'TypeD', 'TypeF_1', 'TypeF_2'};
neuron_IDs = [328533761, 294436967, 792040520, 950229431];
connectivityTable = [];
for iName = 1:numel(neuron_names)
    
    % Identify names of current .csv files
    curr_name = neuron_names{iName};
    us_file_name = [curr_name, '_upstream_partners.csv'];
    ds_file_name = [curr_name, '_downstream_partners.csv'];
    
    % Load upstream partners
    us_partners = readtable(fullfile(parentDir, us_file_name), 'delimiter', ',');
    ID_table = repmat(table({curr_name}, neuron_IDs(iName), {'upstream'}, 'variableNames', ...
            {'DAN_name', 'DAN_ID', 'partnerType'}), size(us_partners, 1), 1);
    us_partners = [ID_table, us_partners];
        
    % Load downstream partners
    ds_partners = readtable(fullfile(parentDir, ds_file_name), 'delimiter', ',');
    ID_table = repmat(table({curr_name}, neuron_IDs(iName), {'downstream'}, 'variableNames', ...
            {'DAN_name', 'DAN_ID', 'partnerType'}), size(ds_partners, 1), 1);
    ds_partners = [ID_table, ds_partners];
    
    % Add both to main table
    connectivityTable = [connectivityTable; us_partners];
    connectivityTable = [connectivityTable; ds_partners];
end
    
save_file = fullfile(parentDir, 'full_connectivity_table.csv');
% writetable(connectivity_table, save_file);

%% LOAD SAVED CONNECTIVTY DATA
parentDir = 'C:\Users\Wilson Lab\Google Drive\Lab Work\LH-DAN_connectomics_analysis';
connectivityTable = readtable(fullfile(parentDir, 'full_connectivity_table.csv'), ...
        'delimiter', ',');
    
%%
dt = DataTable(connectivityTable);
filterDefs = struct();

filterDefs.DAN_name = 'TypeF_2';
filterDefs.partnerType = 'upstream';
filterDefs.type = '';

thresh = 5;

% Get data subset
dt = dt.initialize_filters(filterDefs);
sub = dt.subset;

% Drop rows for partners that are incompletely traced or that do not have an assigned cell type
sub = sub(strcmpi(sub.status, 'traced') & ~cellfun(@isempty, sub.type), :);

% Calculate connection weight of with each partner neuron as a percentage of total DAN synapses
totalSynapses = sum(sub.connectionSynapseCount);
connectionWeightPct = 100 * (sub.connectionSynapseCount ./ totalSynapses);
sub.connectionWeightPct = connectionWeightPct;

% Calculate partner connection weight (percentage of all partner's synapses)
if strcmp(filterDefs.partnerType, 'upstream')
    partnerTotalSynapses = sub.partnerOutputSynapseCount;
elseif strcmp(filterDefs.partnerType, 'downstream')
    partnerTotalSynapses = sub.partnerInputSynapseCount;
end
sub.partnerConnectionWeightPct = 100 * (sub.connectionSynapseCount ./ partnerTotalSynapses);

% Create another table that groups partners by cell type
groupedSub = groupsummary(sub(:, [1 3 5 8:10]), 'type', 'sum', {'connectionSynapseCount', ...
        'partnerInputSynapseCount', 'partnerOutputSynapseCount'});
groupedSub.connectionWeightPct = 100 * (groupedSub.sum_connectionSynapseCount ./ totalSynapses);
if strcmp(filterDefs.partnerType, 'upstream')
    cellTypeTotalSynapses = groupedSub.sum_partnerOutputSynapseCount;
elseif strcmp(filterDefs.partnerType, 'downstream')
    cellTypeTotalSynapses = groupedSub.sum_partnerInputSynapseCount;
end
groupedSub.partnerConnectionWeightPct = 100 * (groupedSub.sum_connectionSynapseCount ./ ...
        cellTypeTotalSynapses);
    
% Plot figures
figure(1);clf; hold on;
plot(sort(groupedSub.connectionWeightPct), '-o');
plot(sort(sub.connectionWeightPct), '-o');
legend({'Cell types', 'Single neurons'}, 'fontsize', 14, 'location', 'nw')
title('Connection weights (% of total LH-DAN input/outputs)')

disp(' '); disp(' ')
disp(tail(sortrows(groupedSub, 3), 10))
disp(tail(sortrows(sub(:, [1 5 8 12]), 4), 10))
disp(tail(sortrows(groupedSub, 7), 10))

figure(2);clf; hold on;
plot(sort(groupedSub.sum_connectionSynapseCount), '-o');
plot(sort(sub.connectionSynapseCount), '-o');
legend({'Cell types', 'Single neurons'}, 'fontsize', 14, 'location', 'nw')
title('Connection synapse counts')

figure(3);clf; hold on;
plot(sort(groupedSub.partnerConnectionWeightPct), '-o');
plot(sort(sub.partnerConnectionWeightPct), '-o')
legend({'Cell types', 'Single neurons'}, 'fontsize', 14, 'location', 'nw')
title('Connection weights (% of total partner neuron input/outputs)')

% Display percentage of partners with connection weight above a certain threshold
nStrongPartners = sum(groupedSub.sum_connectionSynapseCount > thresh);
pctStrongPartners = 100 * (nStrongPartners / numel(groupedSub.sum_connectionSynapseCount));
disp([num2str(round(pctStrongPartners, 1)), '% of ', filterDefs.partnerType, ...
        ' cell types connected by >', num2str(thresh), ' synapses']);
nStrongPartners = sum(sub.connectionSynapseCount > thresh);
pctStrongPartners = 100 * (nStrongPartners / numel(sub.connectionSynapseCount));
disp([num2str(round(pctStrongPartners, 1)), '% of ', filterDefs.partnerType, ...
        ' neurons connected by >', num2str(thresh), ' synapses']);


figure(4);clf; plot(groupedSub.sum_partnerOutputSynapseCount, groupedSub.sum_connectionSynapseCount, '.')
xlabel('partner synapses')
ylabel('connection strength (synapses)')




























    

