function exp_ratio=import_exp_data(filename)
    delimiter = '\t';

    formatSpec = '%f%f%f%f%[^\n\r]';

    fileID = fopen(filename,'r');

    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);

    fclose(fileID);
    exp_ratio = dataArray{:, 3};

    clearvars filename delimiter formatSpec fileID dataArray ans;
end