
testObj = EventAlignedData(load_expList());
testObj = testObj.set_align_event('locomotion');
testObj = testObj.extract_roi_data();
testObj = testObj.extract_fictrac_data();

testObj = testObj.create_filter_event_vectors();
test2 = testObj.output_analysis_subset([2 3]);

% head(test.sourceMd, 1)
% head(test.alignEventTable, 1)
% head(test.eventFlTable, 1)
% head(test.eventFilterTable, 1)

%%

test = test2(strcmp(test2.roiName, 'TypeF-R'), :);
expList = unique(test.expID);
for iExp = 1:numel(expList)
   currExpID = expList{iExp};
   currData = test(strcmp(test.expID, currExpID), :);
   
   currEventFl = cell2mat(currData.eventFl');
   currEventMoveSpeed = cell2mat(currData.moveSpeed');
   
   figure(iExp);clf;
%    plot(smoothdata(currEventFl, 1, 'gaussian', 3));
%    hold on;
%    plot(mean(smoothdata(currEventFl, 1, 'gaussian', 3), 2), 'linewidth', 3, 'color', 'k');
% imagesc(currEventFl');
imagesc(currEventMoveSpeed')
   title(currExpID);
   
end


%%

filts = struct();
filts.roiName = 'TypeD';
filts.ballstop = 0;
filts.flailing = 0;
% filts.grooming = 0;
% filts.isolatedmovement = 0;
% filts.locomotion = 0;
filts.optostim = 0;
filts.odor = 0;
filts.soundstim = 0;
% filts.odorName = 'EtOH';
% filts.concentration = '4%'

filtNames = fieldnames(filts);
filterVec = ones(size(test3, 1), 1);
for iFilt = 1:size(filtNames)
    currFiltName = filtNames{iFilt};
    if isa(filts.(currFiltName), 'char')
        currFiltVec = ~cellfun(@isempty, regexp(test3.(currFiltName), filts.(currFiltName), ...
                'once'));
    elseif filts.(currFiltName) == 0
        currFiltVec = ~cellfun(@any, test3.(currFiltName));
    end
    filterVec = filterVec .* currFiltVec;    
end
filterVec = logical(filterVec);
test4 = test3(filterVec, :);

%

test4 = innerjoin(test4, oneNoteMd(:, {'expID', 'olfactometerVersion'}));

smWin = 3;

expIDList = unique(test4.expID);
subplotDims = numSubplots(numel(expIDList));
if all(subplotDims == [1 3])
   subplotDims = [2 2]; 
end
f = figure(1); clf;
f.Color = [1 1 1];
currPlotNum = 0;
for iExp = 1:numel(expIDList)
    currExpID = expIDList{iExp};
    currData = test4(strcmp(test4.expID, currExpID), :);
    firstRoiName = currData.roiName{1};
    currData = currData(strcmp(currData.roiName, firstRoiName), :);
    currEventFl = cell2mat(currData.eventFl');
    
    fl = currEventFl;
    bl = repmat(currData.trialBaseline', size(fl, 1), 1);
    currEventFl = (fl - bl) ./ bl; 
    
    currEventFl = currEventFl(:, ~any(isnan(currEventFl), 1));
    
    %    currEventMoveSpeed = cell2mat(currData.moveSpeed');
    if size(currEventFl, 2) > 0
        
        currPlotNum = currPlotNum + 1;
        ax = subaxis(subplotDims(1), subplotDims(2), currPlotNum, 'm', 0.02, 'sv', 0.04, 'sh', 0.02);
        
%         plot(smoothdata(currEventFl, 1, 'gaussian', smWin));
        hold on;
        plot(mean(smoothdata(currEventFl, 1, 'gaussian', smWin), 2, 'omitnan'), ...
                'linewidth', 3, 'color', 'k');
%             
        [~, onsetVol] = min(abs(currData.volTimes{1}));  
        stimDur = currData.offsetTime(1) - currData.onsetTime(1);
        [~, offsetVol] = min(abs(currData.volTimes{1} - stimDur));
%         
        
        
        meanData = mean(smoothdata(currEventFl, 1, 'gaussian', smWin), 2, 'omitnan');
        SE = std_err(currEventFl, 2);
        jbfill((1:numel(SE)), (meanData + SE)', (meanData - SE)', [0 0 0], [0 0 0], 1, 0.2);
        yL = ylim();
        
%         plot_stim_shading([onsetVol, offsetVol]);
        plot([onsetVol onsetVol], yL, 'color', 'r', 'linewidth', 3);
        
% imagesc(isnan(currEventFl'));
        ylim(yL);
%                imagesc(currEventFl');
        % imagesc(currEventMoveSpeed')
        title([currExpID, ' - Ver', num2str(currData.olfactometerVersion(1))]);
        box off
    end
end

% unique(test3(:, {'odorName', 'concentration'}))


%%
test = test3(strcmp(test3.expID, '20190606-1'), :);
test = test(contains(test.roiName, 'TypeD-L'), :);
locTest = cell2mat(test.locomotion');
moveSpeedTest = cell2mat(test.moveSpeed');
flTest = cell2mat(test.eventFl');

flTest = flTest(:, ~any(locTest, 1));
figure; imagesc(flTest')
figure; imagesc(moveSpeedTest')
figure; imagesc(locTest');


