function arrTest(arr)
% Displays some basic information about an array (size, min/max, etc).

disp(' ');
disp(['Size: [', num2str(size(arr)), ']']);
disp(['Min: ', num2str(min(arr(:)))]);
disp(['Max: ', num2str(max(arr(:)))]);
if isnan(max(arr(:)))
    disp(['NaN count: ', num2str(sum(isnan(arr(:)))), ' of ' num2str(numel(arr))]);
elseif isinf(max(arr(:)))
    disp(['Inf count: ', num2str(sum(isinf(arr(:)))), ' of ' num2str(numel(arr))]);
end

end