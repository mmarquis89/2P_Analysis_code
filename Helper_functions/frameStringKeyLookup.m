function [ value ] = frameStringKeyLookup( frameString, key )
% From Sasha

cc = strsplit(frameString, '\n');

value = '';
for i=1:length(cc)
    str_tokenized = strtrim(strsplit(cc{i},'='));
    
    if(~isempty(findstr(str_tokenized{1}, key)))
        value = str2num(str_tokenized{2});
    end
end

end

