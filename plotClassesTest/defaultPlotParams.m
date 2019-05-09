classdef defaultPlotParams < plotParams
    
    methods       
        % Constructor
        function obj = defaultPlotParams(goodTrials, varargin)
           if nargin < 2
                groupNames = {'AllTrials'};
           else
               groupNames = varargin{1};
           end
           obj@plotParams(goodTrials);
           obj.groupType = 'Default';
           obj.groupNames = groupNames;
        end
    end
    
    
    
end