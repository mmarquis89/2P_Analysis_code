function cm = custom_colormap(nColors)
% Returns a colormap consisting of nColors visually distinct colors (up to 20, will return entire 
% colormap if called without any arguments)
cm = [rgb('blue'); rgb('red'); rgb('green'); rgb('magenta'); rgb('cyan'); rgb('gold'); ...
    rgb('lime'); rgb('black'); rgb('maroon'); rgb('grey'); rgb('orange'); rgb('olive'); ...
    rgb('darkred'); rgb('indigo'); rgb('darkviolet'); rgb('deepskyblue'); rgb('mediumvioletred'); ...
    rgb('saddlebrown'); rgb('darkgreen'); rgb('salmon')];

if ~isempty(nColors) && nColors <=20
    cm = cm(1:nColors, :);
end

end