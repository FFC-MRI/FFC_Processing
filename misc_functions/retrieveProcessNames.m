function [output] = retrieveProcessNames(as)
%RETRIEVEPROCESSNAMES Summary of this function goes here
%   Detailed explanation goes here

previousLocation = pwd;
cd C:\Users\s02sn2\Desktop\FFCProcessing\FFC-Processing\image_handling\noise_filters
Files = dir('*.m');
num_files = length(Files);
nameArray = {};
nameArray{1} = 'None';


for i=1:num_files
    [pathstr, name, ext] = fileparts(Files(i).name);
    if strcmp(name,'NoiseFilterTemplate') == false
        nameArray{end+1} = feval(name).processName;
    end
end
cd(previousLocation);

output = nameArray;

end

