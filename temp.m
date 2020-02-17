
smWin = 6;

moveCtrlDff = smoothdata(squeeze(bl.dffArr(:, end, :)), 'gaussian', smWin);
moveCtrlRawFl = smoothdata(squeeze(bl.rawFlArr(:, end, :)), 'gaussian', smWin);
moveCtrlZscore = smoothdata(squeeze(bl.zscoreArr(:, end, :)), 'gaussian', smWin);

baseSubFl = moveCtrlRawFl - repmat(min(moveCtrlRawFl), size(moveCtrlRawFl, 1), 1); 

flData = moveCtrlDff;
% flData = moveCtrlRawFl;
% flData = moveCtrlZscore;
flData = baseSubFl;


% plot(sort(baseSubFl(:)));

figure(3);clf;
for iTrial = 1:size(moveCtrlDff, 2)
   subaxis(3, 3, iTrial); hold on;
   currData = flData(:, iTrial);
   currDataBaseSub = currData - movmin(currData, 50);
   plot(currData, 'linewidth', 1);
   plot(currDataBaseSub, 'linewidth', 1);
end

figure(2);clf;
for iTrial = 1:size(moveCtrlDff, 2)
   subaxis(3, 3, iTrial); hold on;
   currData = flData(:, iTrial);
   currDataBaseSub = currData - movmin(currData, 50);
%    histogram(currDataBaseSub, 50);
   plot(sort(currDataBaseSub))
end


% manThresholds = [35 40 40 20 25 25 30 30]; % 2/6 Exp 1
manThresholds = [20 30 20 25 40 30 30 30 30]; % 2/6 Exp 2

goodVols = [];
for iTrial = 1:size(moveCtrlDff, 2)
    currData = flData(:, iTrial);
    currDataBaseSub = currData - movmin(currData, 50);
    goodVols(:, iTrial) = currDataBaseSub < manThresholds(iTrial);
end


