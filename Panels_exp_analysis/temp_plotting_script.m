
expID = '20201124-1';

stimDurs = [60 60 30 60, 60, 60 ,60 ,60];

roiName = 'EB-DAN';

saveFig = 0;

% Adjust plot spacing and margins
SV = 0.05;
SH = 0.05;
ML = 0.05;
MR = 0.07;
MT = 0.05;
MB = 0.08;

f = figure(1);clf; hold on;
f.Color = [1 1 1];
currExpData = roiData(strcmp(roiData.expID, expID) & strcmp(roiData.roiName, roiName), :);

nPlots = max(currExpData.trialNum);

for iTrial = 1:nPlots
    
    ax = subaxis(nPlots, 1, iTrial, 'mt', MT, 'mb', MB, 'sv', SV, 'mr', MR, ...
            'ml', ML, 'sh', SH);
    currTrialData = currExpData(currExpData.trialNum == iTrial, :);
    currTrialData = innerjoin(currTrialData, ftData);
    
    flData = (currTrialData.rawFl{:} - currTrialData.expBaseline) / currTrialData.expBaseline;
    flowData = currTrialData.meanFlow{:};
    
    xx_1 = (1:numel(flData)) ./ 6.87;
    xx_2 = (1:numel(flowData)) .* median(diff(currTrialData.frameTimes{:}));
    
    plot(xx_1, smoothdata(flData, 'gaussian', 5), 'linewidth', 2)
    hold on;
    yL = ylim();
    plot([45, 45], yL, 'color' , 'g', 'linewidth', 4)
    plot([45, 45] + stimDurs(iTrial), yL, 'color' , 'g', 'linewidth', 4)
    ylim(yL)
    ax = gca();
    ax.FontSize = 14;
    if iTrial == nPlots
        xlabel('Time (sec)');
    end
    ylabel('dF/F');
    yyaxis right
    plot(xx_2, smoothdata(flowData, 'gaussian', 50), 'linewidth', 2)
    xlim([0 xx_2(end)])
    ylabel('Optic flow')
    if iTrial == 1
       title(expID) 
    end

end

if saveFig
    figDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\EB-DAN_GroupedAnalysisData\Figs';
    save_figure(f, figDir, [expID, '_single_trial_dff+optic_flow'])
end


