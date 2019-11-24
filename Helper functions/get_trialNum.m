function trialNum = get_trialNum(str)
    % Extracts a trial number from a string or cell array of strings containing padded trial 
    % numbers (padded to 3 digits). Returns either a scalar or a numeric vector with the results.
    % ex. 'myFile_trial_001.mat' --> 1
    % ex. {'myFile_trial_001.mat', 'myFile_trial_002.mat'} --> [1, 2]
    
    trialNum = str2double(regexp(str, '(?<=trial_)...', 'match', 'once'));

end