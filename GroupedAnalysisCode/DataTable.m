classdef DataTable
% ==================================================================================================   
%    
% Properties:
%       sourceData
%       filterDescriptions
%       filters
%       fieldNames (dependent)
%
% Methods:
%       .clear_all_filters()
%       .clear_filter(fieldName)
%       .add_filter(fieldName, filter)
%       .apply_filters() 
%       .list_filters()
%       .show_filters() 
% 
% ==================================================================================================
    properties
        sourceData          % source table
        filterDescriptions  % description of current filter status 
        filters             % logical array of filter vectors, same dimensions as sourceData 
    end
    properties (Dependent)
        fieldNames          % field names of the source data table
    end 
    methods
        
        % Contructor
        function obj = dataTable(sourceData)
           obj.sourceData = sourceData;
           obj.filters = logical(ones(size(sourceData)));
           obj.filterDescriptions = cell(1, numel(fieldnames(obj.sourceData)));
        end
        
        % Clear all filters
        function obj = clear_all_filters(obj)
            obj.filters = logical(ones(size(obj.sourceData)));
            obj.filterDescriptions = cell(1, numel(fieldnames(obj.sourceData)));
        end
        
        % Clear filter for single field
        function obj = clear_filter(obj, fieldName)
            fNames = fieldnames(obj.sourceData);
            fieldIndex = find(strcmp(fNames, fieldName));
            obj.filterDescriptions{fieldIndex} = [];
            obj.filters(:, fieldIndex) = 1;
        end
        
        % Add a new filter
        function obj = add_filter(obj, fieldName, filter)
           
            fNames = fieldnames(obj.sourceData);
            fieldIndex = find(strcmp(fNames, fieldName));
            
            % Interpret filter according to variable type
            if isa(filter, 'function_handle')
                % Custom function to apply to a column of data
                obj.filterDescriptions{fieldIndex} = filter;
                obj.filters(:, fieldIndex) = filter(obj.sourceData.(fieldIndex));
            elseif isa(filter, 'char')
                % String
                obj.filterDescriptions{fieldIndex} = ['equals "', filter, '"'];
                obj.filters(:, fieldIndex) = strcmp(obj.sourceData.(fieldIndex), filter);
            elseif isa(filter, 'logical') 
                % Logical vector
                obj.filterDescriptions{fieldIndex} = 'User-provided logical vector';
                obj.filters(:, fieldIndex) = filter;
            elseif isa(filter, 'numeric') 
                if numel(filter) == size(obj.sourceData, 1) 
                    % Numeric vector with length == table row count
                    obj.filterDescriptions{fieldIndex} = 'User-provided logical vector';
                    obj.filters(:, fieldIndex) = logical(filter);
                else
                    % Scalar or numeric vector with length < table row count
                    obj.filterDescriptions{fieldIndex} = ['equals any of: ', num2str(filter)];
                    obj.filters(:, fieldIndex) = ismember(obj.sourceData.(fieldIndex), filter);
                end
            elseif isa(filter, 'cell')
                % Cell array
                obj.filters(:, fieldIndex) = ismember(obj.sourceData.(fieldIndex), filter);
                filterStr = 'is any of: {';
                for iStr = 1:numel(filter)
                   filterStr = [filterStr, filter{iStr}, ' ']; 
                end
                obj.filterDescriptions{fieldIndex} = [filterStr(1:end-1), '}'];
            end
        end
        
        % Return filtered source table
        function outputData = apply_filters(obj)
            finalFilter = ones(size(obj.sourceData, 1), 1);
            for iField = 1:size(obj.sourceData, 2)
                finalFilter = finalFilter .* obj.filters(:, iField);
            end
            outputData = obj.sourceData(logical(finalFilter), :);
        end
        
        % List filter descriptions for each field
        function list_filters(obj)
                
            dispTable = table(fieldnames(obj.sourceData), obj.filterDescriptions', ...
                    'VariableNames', {'Field_name', 'Filter'});
            disp(dispTable)
        end
        
        % Show 2D plot of current filter status
        function show_filters(obj)
            figure(1);clf; imagesc(obj.filters);
        end
        
        % Return source table field names, formatted in a view-friendly table
        function fieldNames = get.fieldNames(obj)
           rawNames = fieldnames(obj.sourceData);
           col = (1:numel(rawNames))';
           fieldNames = (table(col, rawNames)); 
        end
        
    end



end