classdef moveSpeedFilteredPlotParams < plotParams
    
    properties
        sdCap
        frameRate
        smWin
        smSpeedData
    end
    
    methods
        % Constructor
        function obj = moveSpeedFilteredPlotParams(goodTrials, ftData, varargin)
            
            obj@plotParams(goodTrials);
            
            p = inputParser;
            addParameter(p, 'sdCap', 5);
            addParameter(p, 'FrameRate', 25);
            addParameter(p, 'SmoothWinSize', 5);
            parse(p, varargin{:});
            obj.sdCap = p.Results.sdCap;
            obj.frameRate = p.Results.FrameRate;
            obj.smWin = p.Results.SmoothWinSize;
            
            obj.groupType = '';
            obj.groupNames = [];
            obj.groupShading = [];
            
            % Extract FicTrac speed data
            rawData = ftData.moveSpeed;                    % --> [frame, trial]
            rawData = rawData';                            % --> [trial, frame]
            plotData = (rad2deg(rawData .* obj.frameRate));   % --> [trial, frame] (deg/sec)
            
            % Cap values at n SD above mean
            capVal = mean(plotData(:), 'omitnan') + (obj.sdCap * std(plotData(:), 'omitnan'));
            plotData(plotData > capVal) = capVal;
            
            % Smooth data
            obj.smSpeedData = smoothdata(plotData, 2, 'gaussian', obj.smWin);
            
        end
        
        % Create high vs. low speed trial groups
        function obj = make_highVsLowSpeedGroups(obj, )
            
            
        end
        
        
        
    end
    
    
    
end