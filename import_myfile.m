function outgenenamefile_list=import_myfile(filename)

    delimiter = '';

    formatSpec = '%s%[^\n\r]';

    fileID = fopen(filename,'r');

    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);

    fclose(fileID);

    outgenenamefile_list = dataArray{:, 1};


    clearvars filename delimiter formatSpec fileID dataArray ans;