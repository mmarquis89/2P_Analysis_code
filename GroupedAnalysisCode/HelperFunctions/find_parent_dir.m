function parentDir = find_parent_dir(expID)
% Determine parent imaging data dir based on an experiment date

expDateNum = str2double(expID(1:8));
if expDateNum < 20180206
    parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017';
elseif expDateNum < 20180405
    parentDir = 'F:\ImagingData';
elseif expDateNum < 20180623
    parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2018\2018 Apr-May';
elseif expDateNum < 20190211
    parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2018';
else
    error('Not ready to process those experiments yet');
end

end