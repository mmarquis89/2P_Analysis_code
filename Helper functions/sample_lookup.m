classdef sample_lookup
% Converts between sample indices of data with different sampling rates
% e.g. 
%   vol2frame = sample_lookup(volumeRate, frameRate);
%   stimStartFrames = frame2vol.convert(stimStartVols);  
% 
    properties
        baseSampRate
        querySampRate
    end
    
    methods
        
        % Constructor
        function obj = sample_lookup(convertToSampRate, convertFromSampRate)
            obj.baseSampRate = convertToSampRate;
            obj.querySampRate = convertFromSampRate;            
        end
        
        % Sample lookup/conversion
        function baseIdxMatch = convert(obj, queryIdx)
            querySampleTimes = queryIdx / obj.querySampRate;
            maxBaseSample = ceil(max(querySampleTimes) * obj.baseSampRate);
            baseSamples = 1:maxBaseSample;
            baseSampleTimes = baseSamples / obj.baseSampRate;
            baseIdxMatch = interp1(baseSampleTimes, baseSamples, querySampleTimes, 'nearest');
            if querySampleTimes(1) < baseSampleTimes(1)
                baseIdxMatch(1) = 1;
            end
            if querySampleTimes(end) > baseSampleTimes(end)
               baseIdxMatch(end) = maxBaseSample; 
            end
        end
    end
end