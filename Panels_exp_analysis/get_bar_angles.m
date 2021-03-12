function barAngles = get_bar_angles(Pats)
%===================================================================================================
% Returns a vector containing the bar angle in degrees relative to the front of the fly (+- 180, 
% positive == left) for each X-index of the panels pattern array. This vector can then be used to 
% convert any bar position data from pattern indices to angles in degrees.
%===================================================================================================

% Drop Y-dimension of pattern array if it exists
if numel(size(Pats)) > 3
    Pats = Pats(:, :, :, 2);
end

% Get angles in degrees of each LED pixel on the panels (arranged by pattern index)
pxSize = 3.75;
startAngle = -180 + (pxSize * 4.5);
panelsPixelAngles = [startAngle:pxSize:0, ...
                    (pxSize/2):pxSize:180, ...
                    -180 + (pxSize/2):pxSize:(startAngle - pxSize)];
                
% Locate the midpoint of the bar in each pattern index
barAngles = nan(size(Pats, 3), 1);
for iPos = 1:size(Pats, 3)
    leftEdge = find(Pats(end, :, iPos), 1, 'first');
    rightEdge = find(Pats(end, :, iPos), 1, 'last'); 
    if ~isempty(leftEdge)
        if  abs(panelsPixelAngles(leftEdge)) > 90 ...
                && panelsPixelAngles(leftEdge) == -panelsPixelAngles(rightEdge)
            % Prevent error if bar is directly behind the fly
            barAngles(iPos) = 180;
        else
            barAngles(iPos) = mean(panelsPixelAngles([leftEdge, rightEdge]));
        end
    end
end  

end