function [vectorStrength, vectorPhase] = calculate_vector_strength(phaseData, respData)
%===================================================================================================
% Calculate the vector strength and phase from a continuous response (e.g. voltage or GCaMP 
% fluorescence) to a circular variable (e.g. visual stim position on panels).
%
%  Analysis modified from YF and originally based on:?
%       Xiang Gao & Wehr, M. Neuron 86, 292–303 (2015)
%       ?Zar, JH. Biostatistical analysis. 4. Upper Saddle River, N.J: Prentice Hall; 1999.
%
% INPUTS:
%
%       phaseData   = vector of phase data with values ranging from 0 to 2*pi
% 
%       respData    = the response data corresponding to each point in the phaseData vector (must 
%                     be same dimensions as phaseData)
% 
% OUTPUTS:
%
%       vectorStrength  = the calculated vector strength in respData units of magnitude
% 
%       vectorPhase     = the calculated vector phase in terms of the input phase data
% 
% 
%===================================================================================================

n = numel(phaseData);

opposite = sum(respData .* sin(phaseData));
adjacent = sum(respData .* cos(phaseData));
oppositeMean =  opposite / n;
adjacentMean =  adjacent / n;

vectorStrength = sqrt(opposite^2 + adjacent^2) / n;
vectorPhase = wrapTo2Pi(atan2(oppositeMean, adjacentMean)); % atan2() vs atan() is important here


end