function genename = import_gene_name(filename, startRow, endRow)

%
% Example:
%   VarName0 = importfile1('0.1_MF_ATP_gene_name.txt',1, 3959);
%

delimiter = '';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

formatSpec = '%s%[^\n\r]';

fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end

fclose(fileID);


genename = dataArray{:, 1};


