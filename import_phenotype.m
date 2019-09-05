function phenotype = import_phenotype(filename, startRow, endRow)

% Example:
%   phenotype1 = importfile1('phenotype.txt', 1, 3);
%

delimiter = '\t';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';

fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

fclose(fileID);


phenotype = [dataArray{1:end-1}];



