classdef DataTable
% ==================================================================================================   
%    
% Properties:
%       sourceData
%       filterDefs
%       filterMat
%       filterVec
%       fieldNames (dependent)
%       subset (dependent)
% Methods:
%       .generate_filters(filterDefs)
%       .update_filter(filtName, filtValue)
%       .show_filters() 
% 
% ==================================================================================================
properties
    sourceData          % source table
    filterDefs          % struct containing filters for fields of sourceData
    filterMat           % logical array of filter vectors, same dimensions as sourceData
    filterVec           % Single logical vector for table, combining all filters
    customProps         % Additional properties of the source data
end

properties (Dependent)
    fieldNames          % field names of the source data table
    subset              % a version of sourceData, using filterVec to subset rows
end

methods
    
    % Contructor
    function obj = DataTable(sourceTable)
        obj.sourceData = sourceTable;
        obj.filterMat = true(size(sourceTable));
        obj.filterVec = true(size(sourceTable, 1), 1);
        if ~isempty(sourceTable.Properties.CustomProperties) %#ok<*MCNPN>
            fNames = fieldnames(sourceTable.Properties.CustomProperties);
            for iField = 1:numel(fNames)
                currField = fNames{iField};
                obj.customProps.(currField) = sourceTable.Properties.CustomProperties.(currField);
            end
        end
    end
    
    % Add a new filter
    function obj = generate_filters(obj, filterDefs)
%   if filter is a function handle: 
%       - Apply the function to any field and return the result (must be n x 1 logical vector)    
%
%   If filter is a string:
%       - Apply a regex string to a character field
%           Example: 'regex'
% 
%       - Evaluate as an expression to apply to a numeric field, replacing 'x' with 'filtData'
%           Example: '(x > 4 | x == 2) & x ~= 20'
% 
%   If filter is numeric:
%       - Select a specific value or set of values from a numeric scalar field:
%           Example: 2 or [2 4 6]
% 
%       - Specify a code for a logical condition to apply to a logical event vector field:
%           0:  remove events with any occurance of the filter event
%           -1: remove events with any pre-onset occurance of the filter event
%           -2: remove events with any post-onset occurance of the filter event
%           1: remove events withOUT any pre-onset occurance of the filter event
%           2: remove events withOUT any post-onset occurance of the filter event
        
    % Reset filterMat and filterVec
    obj.filterMat = true(size(obj.sourceData));
    obj.filterVec = true(size(obj.sourceData, 1), 1);
    
    % Create new filters
        filtNames = fieldnames(filterDefs);

        for iFilt = 1:numel(filtNames)
            
            % Get current filter info
            filtName = filtNames{iFilt};
            filtValue = filterDefs.(filtName);
            filtData = obj.sourceData.(filtName);
            
            % ------ Create the filter ------
            if ~isempty(filtValue)
                
                currFiltVec = parse_filter(obj, filtValue, filtData);
                
                % Update appropriate column of filterMat
                obj.filterMat(:, strcmp(fieldnames(obj.sourceData), filtName)) = currFiltVec;

            end
        end%iFilt
        
        % Update master filterVec
        obj.filterVec = all(obj.filterMat, 2);
        
        obj.filterDefs = filterDefs;
    end
    
    % Update a single filter
    function obj = update_filter(obj, filtName, filtValue)
        
        % Create the new filter vector
        filtData = obj.sourceData.(filtName);        
        newFiltVec = parse_filter(obj, filtValue, filtData);
        
        % Update the appropriate column in the filter matrix
        obj.filterMat(:, strcmp(fieldnames(obj.sourceData), filtName)) = newFiltVec;
        
        % Update the main filterVec
        obj.filterVec = all(obj.filterMat, 2);
        
        % Update filterDefs
        obj.filterDefs.(filtName) = filtValue;
        
    end
    
    % Show 2D plot of current filter status
    function show_filters(obj)
        figure(1); clf;
        ax = gca();
        imagesc(obj.filterMat);
        ax.XTick = 1:size(obj.filterMat, 2);
        ax.XTickLabel = obj.sourceData.Properties.VariableNames;
        ax.XTickLabelRotation = 60;
        ax.XAxisLocation = 'top';
        
    end
    
    % Return filtered source table
    function outputTable = get.subset(obj)
        
        % Create output table
        outputTable = obj.sourceData(logical(obj.filterVec), :);
        
        % Copy any custom props over to the output table
        if ~isempty(obj.customProps)
            fNames = fieldnames(obj.customProps);
            outputTable = addprop(outputTable, fNames, repmat({'table'}, size(fNames)));
            for iField = 1:numel(fNames)
                currField = fNames{iField};
                outputTable.Properties.CustomProperties.(currField) = obj.customProps.(currField);
            end
        end
    end
    
    % Return source table field names, formatted in a view-friendly table
    function fieldNames = get.fieldNames(obj)
        rawNames = obj.sourceData.Properties.VariableNames;
        col = (1:numel(rawNames))';
        fieldNames = (table(col, rawNames'));
    end
    
end%methods

end%class

% ==================================================================================================
% Local functions
% ==================================================================================================

% Determine data type of a particular filterDef field value
function filterDataType = get_filter_type(filtValue)
if isa(filtValue, 'char')
    filterDataType = 'charVector';
elseif isa(filtValue, 'function_handle')
    filterDataType = 'functionHandle';
elseif isa(filtValue, 'numeric') && isscalar(filtValue)
    filterDataType = 'numericScalar';
else
    filterDataType = 'numericVector';
end
end

% Determine data type of a particular table field to use in filtering
function fieldDataType = get_data_type(filtData)
i = 1;
while isempty(filtData{i})
    i = i + 1;
end
if isa(filtData, 'cell') && ischar(filtData{i})
    fieldDataType = 'charVector';
elseif isa(filtData, 'cell') && isnumeric(filtData{i})
    fieldDataType = 'numericVector';
else
    fieldDataType = 'numericScalar';
end
end

% Decide how to interpret a filter and return the appropriate logical vector
function currFiltVec = parse_filter(obj, filtValue, filtData)

% Check for alignment event name
if ismember(fieldnames(obj.customProps), 'alignEventName')
    alignEventName = obj.customProps.alignEventName;
else
    alignEventName = [];
end

% Determine the data type of the filter value and the table field to be filtered
filterType = get_filter_type(filtValue);
dataType = get_data_type(filtData);

% Custom function to apply to a column of data
if strcmp(filterType, 'functionHandle')
    currFiltVec = filtValue(filtData);
    
    % A regex string for a char field
elseif strcmp(filterType, 'charVector') && strcmp(dataType, 'charVector')
    test = regexp(filtData, filtValue, 'once');
    currFiltVec = ~cellfun(@isempty, test);
%     currFiltVec = ~cellfun(@isempty, regexp(filtData, filtValue, 'once'));
    
    % A string to be evaluated as an expression after replacing 'x' with 'filtData'
elseif strcmp(filterType, 'charVector')
    currFiltVec = eval(regexprep(filtValue, 'x', 'filtData'));
    
    % A value specifying a logical filter vector condition
elseif strcmp(filterType, 'numericScalar') && strcmp(dataType, 'numericVector')
    
    onsetVols = cellfun(@(x) argmin(abs(x)), obj.sourceData.volTimes);
    onsetVols = mat2cell(onsetVols, ones(size(onsetVols)));
    
    if ismember(fieldnames(obj.customProps), 'alignEventName') & strcmp(filtName, ...
            alignEventName)
        
        % Filter out events with a pre-onset occurance of the alignment event type
        if filtValue < 1
            
            currFiltVec = ~cellfun(@(x, y) any(x(1:(y - 1))), filtData, onsetVols);
            
            % Filter out events withOUT a pre-onset occurance of the alignment event time
        elseif filtValue == 1
            
            currFiltVec = cellfun(@(x, y) any(x(1:(y - 1))), filtData, onsetVols);
        end
        
        % Filter out events with any occurance of the filter event type
    elseif filtValue == 0
        currFiltVec = ~cellfun(@any, filtData);
        
        % Filter out events with any pre-onset occurance of the filter event type
    elseif filtValue == -1
        currFiltVec = ~cellfun(@(x, y) any(x(1:(y - 1))), filtData, onsetVols);
        
        % Filter out events withOUT any pre-onset occurance of the filter event type
    elseif filtValue == 1
        currFiltVec = cellfun(@(x, y) any(x(1:(y - 1))), filtData, onsetVols);
        
        % Filter out events with any post-onset occurance of the filter event type
    elseif filtValue == -2
        currFiltVec = ~cellfun(@(x, y) any(x(y:end)), filtData, onsetVols);
        
        % Filter out events withOUT any post-onset occurance of the filter event type
    elseif filtValue == 2
        currFiltVec = cellfun(@(x, y) any(x((y:end))), filtData, onsetVols);
        
    else
        error(['Invalid filter type "', filterType, '"']);
    end
    
    % A specific scalar value or set of acceptable scalar values
elseif strcmp(dataType, 'numericScalar')
    currFiltVec = ismember(filtData, filtValue);
    
elseif isempty(filtValue)
    currFiltVec = ones(size(filtData));
else
    error(['Invalid filter type "', filterType, '"']);
end
end
