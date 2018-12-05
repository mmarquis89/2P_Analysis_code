function peek(dataArr)
%===================================================================================================
% 
% Creates or activates figure #10, clears the screen and plots the 2D dataArr using imagesc. This 
% function's only purpose is to let me do this without much typing.
%
% INPUTS:
%       onsetArr  = 2D array of numeric data that you want to visualize.
%
%===================================================================================================

figure(10);
clf;
imagesc(dataArr);

end