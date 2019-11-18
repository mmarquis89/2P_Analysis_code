function trialNum = get_trialNum(str)
    % Extracts a trial number from a string containing a padded trial number.
    % ex. "myFile_trial_001.mat" --> 1
    
    trialNum = str2double(regexp(str, '(?<=trial_)...', 'match', 'once'));

end