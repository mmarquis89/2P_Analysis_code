function cm = custom_colormap(nColors)
% Returns a colormap consisting of nColors visually distinct colors (up to 20, will start over again
% at the beginning if nColors > 20). Returns entire colormap if called without any arguments.
cm = [rgb('blue'); rgb('red'); rgb('green'); rgb('magenta'); rgb('cyan').* 0.95; rgb('gold'); ...
    rgb('lime').* 0.95; rgb('black'); rgb('maroon'); rgb('grey'); rgb('orange'); rgb('olive'); ...
    rgb('darkred'); rgb('indigo'); rgb('darkviolet'); rgb('deepskyblue'); rgb('mediumvioletred'); ...
    rgb('saddlebrown'); rgb('darkgreen'); rgb('salmon')];

if ~isempty(nColors) 
    if nColors > 20
        cm = repmat(cm, ceil(nColors / 20), 1);
    end
    cm = cm(1:nColors, :);
end

end