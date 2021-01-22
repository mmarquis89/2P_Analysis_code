function volTime = sec2vols(secTime, volumeRate)

% Converts a time from units of seconds to imaging volumes
volTime = floor(secTime * volumeRate);

end