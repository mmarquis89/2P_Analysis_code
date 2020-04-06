function meanFlowMag = optic_flow_calc(vidFile, varargin)
% %==================================================================================================
% CALCULATE MEAN OPTIC FLOW DATA
%
% Calculates the mean optic flow magnitude for each frame of a video file. Can optionally restrict 
% the flow calculation to just a particular ROI by passing a mask with the same dimensions as the 
% input video frames. NOTE: the mean flow value for the first frame is undefined and so is set to 
% zero.
%
% INPUTS:
%
%       vidFile     = the full path of a video file to use for flow calculation
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%
%       'roiMask'   = a binary array with the same dimensions as the input video frames specifying 
%                     an ROI to focus on. All pixels that are zeros in the mask will be ignored for 
%                     the averaging of the flow data.
%
% OUTPUTS:
%
%       meanFlowMag = an vector with the average flow value for each frame in the video
% 
%===================================================================================================       

% Parse optional arguments
p = inputParser;
addParameter(p, 'roiMask', []);
parse(p, varargin{:});
roiMask = p.Results.roiMask;

% Load video file
vidObj = VideoReader(vidFile);

% Calculate optic flow
opticFlow = opticalFlowFarneback;
frameCount = 0;
meanFlowMag = []; 
while hasFrame(vidObj)
    frameCount = frameCount + 1;
    if ~mod(frameCount, 1000)
       disp(['Processing frame ', num2str(frameCount)]); 
    end
    currFrame = readFrame(vidObj);
    currFrame = currFrame(:,:,1);
    currFlow = estimateFlow(opticFlow, currFrame);
    currFlow = currFlow.Magnitude;
    if ~isempty(roiMask)
        currFlow(~roiMask) = nan;        
    end
    meanFlowMag(end + 1) = mean(currFlow(:), 'omitnan');
end
meanFlowMag(1) = 0;
disp(['Processed ', num2str(frameCount), ' video frames']);
end







