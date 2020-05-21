function outputTable = inner_join(varargin)
% ==================================================================================================   
%  Performs an inner join operation on 3+ tables by calling the "innerjoin" repeatedly on all 
%  arguments.
%
%  INPUTS: 
%       table_1, table_2, ...
%
%  OUTPUTS:
%       outputTable
%
% ==================================================================================================

outputTable = varargin{1};
for iArg = 2:nargin
    outputTable = innerjoin(outputTable, varargin{iArg});    
end

end