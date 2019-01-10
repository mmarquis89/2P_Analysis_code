function remake_block_vid(vidDataDir, blockVidName, sid, bid)
try
   % FlyCap mjpeg AVI files seem to have broken indexing when they exceed a certain file size, so 
   % this function uses ffmpeg to make a (renamed) copy of the file with corrected indexing.
   sourceVid = fullfile(vidDataDir, blockVidName);
   outputVidName = ['block_vid_sid_', num2str(sid), '_bid_', pad(num2str(bid), 3, 'left', '0')];
   outputVid = fullfile(vidDataDir, [outputVidName, '.avi']);
   cmdStr = ['module load ffmpeg/3.3.3; ffmpeg -i "', sourceVid, '" -codec copy "', outputVid, '"'];
   
   % Run only if these output files don't already exist in the output directory
   if ~exist(outputVid, 'file')
       system(cmdStr, '-echo');
   else
       write_to_log(['Output file ', outputVid, ' already exists...skipping vid remake'], mfilename);
   end

catch ME
    write_to_log(getReport(ME), mfilename);
end