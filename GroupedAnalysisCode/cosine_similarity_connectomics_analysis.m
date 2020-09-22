%% LOAD SAVED CONNECTIVTY DATA

parentDir = 'C:\Users\Wilson Lab\Google Drive\Lab Work\LH-DAN_connectomics_analysis';
tbl = readtable(fullfile(parentDir, 'cosine_similarity_analysis_data_supplement_all_MB_DANs.csv'), ...
        'delimiter', ',');

% Convert neuronIDs column from numeric to string
neuronIDs = {};
for i=1:size(tbl, 1)
    disp(i)
    neuronIDs{i} = num2str(tbl.neuronID(i));
end
tbl.neuronID = neuronIDs';

%%

tbl_MBON = tbl(contains(tbl.neuronType, 'MBON'), :);

tbl_1 = readtable(fullfile(parentDir, 'cosine_similarity_analysis_data_all_DANs.csv'), 'delimiter', ',');
% tbl_2 = readtable(fullfile(parentDir, 'cosine_similarity_analysis_data_supplement_AL-PNs.csv'), 'delimiter', ',');
% tbl_3 = readtable(fullfile(parentDir, 'cosine_similarity_analysis_data_supplement_LHNs.csv'), 'delimiter', ',');
% tbl_4 = readtable(fullfile(parentDir, 'cosine_similarity_analysis_data_supplement_PPL1-2+PPM12.csv'), 'delimiter', ',');

tbl = [tbl_MBON; tbl_1];
tbl = unique(tbl);

%%

% Generate unique names for each neuron and convert to a table
newList = {};
count = 1;
currName = '';
for iNeuron = 1:size(neuronList, 1)
    oldName = currName;
    currName = neuronList.neuronName{iNeuron};
    if strcmp(oldName, currName)
        count = count + 1;
    else
        count = 1;
    end
    newList{iNeuron, 1} = [currName, '-', num2str(count)];
end

% Convert CsMat into a table with named columns and rows
csTbl = table('size', size(CsMat), 'VariableNames', newList, 'RowNames', newList, ...
        'VariableTypes', repmat({'double'}, size(CsMat, 1), 1));
for iCol = 1:size(csTbl, 2)
    csTbl{:, iCol} = CsMat(:, iCol);
end

% save(fullfile(parentDir, 'csTbl.mat'), 'csTbl');

%% Initialize data
dt = DataTable(tbl);
filterDefs = struct();
filterDefs.partnerDirection = 'downstream';
filterDefs.status = 'Traced';
filterDefs.partnerCropped = 'FALSE';
dt = dt.initialize_filters(filterDefs);

% Get list of all query neurons
neuronList = unique(dt.subset(:, {'neuronID', 'neuronName'}));
neuronList = sortrows(neuronList, 2);

% Calculate pairwise cosine similarity
CsMat = [];
for iType = 1:size(neuronList, 1)
    currNeuron = neuronList.neuronID{iType};
    disp([num2str(iType), ' of ', num2str(size(neuronList, 1))])
    for iCompType = 1:size(neuronList, 1)
%         if iType ~= iCompType
            currCompNeuron = neuronList.neuronID{iCompType};
            CsMat(iType, iCompType) = compute_similarity(dt.subset, currNeuron, currCompNeuron);
%         end
    end
end


save(fullfile(parentDir, 'CsMat_allDANs+MBONs_downstream.mat'), 'CsMat', 'neuronList');

filterDefs.partnerDirection = 'upstream';
filterDefs.status = 'Traced';
filterDefs.partnerCropped = 'FALSE';
dt = dt.initialize_filters(filterDefs);

% Get list of all query neurons
neuronList = unique(dt.subset(:, {'neuronID', 'neuronName'}));
neuronList = sortrows(neuronList, 2);

% Calculate pairwise cosine similarity
CsMat = [];
for iType = 1:size(neuronList, 1)
    currNeuron = neuronList.neuronID{iType};
    disp(iType)
    for iCompType = 1:size(neuronList, 1)
%         if iType ~= iCompType
            currCompNeuron = neuronList.neuronID{iCompType};
            CsMat(iType, iCompType) = compute_similarity(dt.subset, currNeuron, currCompNeuron);
%         end
    end
end


save(fullfile(parentDir, 'CsMat_allDANs+MBONs_upstream.mat'), 'CsMat', 'neuronList');

% 
%%

allNames = csTbl.Properties.RowNames;

currNeuronInds = contains(csTbl.Properties.RowNames, 'PPL2') | ...
                contains(csTbl.Properties.RowNames, 'PPM12') | ...
                contains(csTbl.Properties.RowNames, 'PAL') | ...
                contains(csTbl.Properties.RowNames, 'PPL1') | ...
                contains(csTbl.Properties.RowNames, 'PAM');
currNeuronInds(contains(csTbl.Properties.RowNames, 'L-1')) = 0;

            
test = csTbl(currNeuronInds, currNeuronInds);

figure(2);clf;
imagesc(test{:, :})
ax = gca;
ax.XTick = 1:size(test, 2);
ax.YTick = 1:size(test, 1);
ax.XTickLabel = test.Properties.VariableNames;
ax.XTickLabelRotation = 90;
ax.YTickLabel = test.Properties.RowNames;
axis square

%%
function Cs = compute_similarity(tbl, neuron_1, neuron_2)

% Separate data for target neurons
sb_1 = tbl(strcmp(tbl.neuronID, neuron_1), :);
sb_2 = tbl(strcmp(tbl.neuronID, neuron_2), :);

% Join synapse count data for all partners
allPartners = outerjoin(sb_1(:, {'partnerID', 'synapseCount'}), sb_2(:, ...
        {'partnerID' 'synapseCount'}), 'Keys', 'partnerID', 'MergeKeys', 1);

% Fill in zero values for non-shared connections
allPartners.synapseCount_left(isnan(allPartners.synapseCount_left)) = 0;
allPartners.synapseCount_right(isnan(allPartners.synapseCount_right)) = 0;

Cs = getCosineSimilarity(allPartners.synapseCount_left, allPartners.synapseCount_right);
    
end