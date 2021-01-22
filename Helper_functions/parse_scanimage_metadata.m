function scanImageInfo = parse_scanimage_metadata(tifMetadata)
%===================================================================================================
% 
% Takes a structure of the metadata extracted from a ScanImage .tif file and reorganizes it to make 
% it easier to find relevant information. This structure should be the metadata output of the 
% read_patterned_tifdata (and therefore also the read_tif) function. Here is a description of
% the reorganization process:
%
% - All the fields from the top level except [tifdata] are combined with all the fields within 
%   [tifdata] except [Software] and added to a field named [fileInfo] in the output structure
%
% - The contents of the [Software] field are parsed and converted to a structure, which is then 
%   added to the output structure as a field called [SI].
%
% 
% ==================================================================================================

mD = tifMetadata; % for conciseness

% To handle files created by both old and new versions of ScanImage
scanImageInfo.fileInfo = rmfield(mD, 'tifinfo'); 
try
    scanImageInfo.fileInfo = setstructfields(scanImageInfo.fileInfo, ...
            rmfield(mD.tifinfo, 'Software'));
catch
    scanImageInfo.fileInfo = setstructfields(scanImageInfo.fileInfo, ...
            rmfield(mD.tifinfo, 'ImageDescription'));
end

% Split up ScanImage data string
try
    siDataStr = mD.tifinfo.Software;
catch
   siDataStr = mD.tifinfo.ImageDescription; 
end

siDataCell = strsplit(siDataStr, '\n');

% Loop through each entry and convert the string to a field in the output structure
scanImageInfo.SI = [];
for iParam = 1:numel(siDataCell)
    currStr = siDataCell{iParam};
    if ~isempty(currStr)
       
       % The entries for each param are formatted as code, so we can just use the eval() function 
       try
%            disp(currStr)
           eval([regexprep(currStr, '.*\.?SI\.', 'scanImageInfo.SI.'), ';'])
       catch
           disp(['Error evaluating "', currStr, '"'])
       end
    end
end

end%function