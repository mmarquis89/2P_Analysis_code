

exePath = 'C:\My Programs\ffmpeg\ffmpeg-20170807-1bef008-win64-static\bin\ffmpeg.exe';

vidDir = 'D:\Dropbox (HMS)\2P Data\Behavior Vids\2018_12_03_exp_1';
inputVidName = 'fc2_save_2018-12-03-151812-0000.avi'
outputVidName 

cmdStr = ['ffmpeg -i "', fullfile(vidDir, inputVidName), '" -codec copy fc2_new.avi']


system(cmdStr, '-echo');






%%
expDates = { ...
    '2018_12_03_exp_2', ...
    '2018_12_03_exp_3', ...
    '2018_12_03_exp_4', ...
    '2018_12_04_exp_1', ...
    '2018_12_04_exp_2', ...
    '2018_12_04_exp_' ...
    };

sid = 0;
FRAME_RATE = 25;

for iExp = 1:numel(expDates)
    
    expDate = expDates{iExp};
    
    parentDir = fullfile('D:\Dropbox (HMS)\2P Data\Behavior Vids', expDate);
    
    % Identify video files
    vidFiles = dir(fullfile(parentDir, '*.mp4'));
    nBlocks = numel(vidFiles);
    for iBlock = 1:nBlocks
        
        % Open VideoReader and VideoWriter for current file
        myVid = VideoReader(fullfile(parentDir, vidFiles(iBlock).name));
        outputVidName = ['sid_', num2str(sid), '_bid_', num2str(iBlock-1), '.avi'];
        myVidWriter = VideoWriter(fullfile(parentDir, outputVidName), 'Motion JPEG AVI');
        myVidWriter.FrameRate = FRAME_RATE;
        open(myVidWriter)
        frameCount = 0;
        while hasFrame(myVid)
            if ~mod(frameCount, 100)
                disp(['Processing frame ', num2str(frameCount), ' of ', ...
                    num2str(FRAME_RATE * myVid.Duration)])
            end
            if ~mod(frameCount, 1000)
                disp(['Working on ', expDate, ', block ', num2str(iBlock)])
            end
            frameCount = frameCount + 1;
            currFrame = readFrame(myVid);
            writeFrame = uint8(currFrame(:,:,1));
            writeVideo(myVidWriter, currFrame)
        end
        close(myVidWriter)
    end%iBlock
end%iExp