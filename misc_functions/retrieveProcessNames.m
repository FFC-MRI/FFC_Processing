function [output] = retrieveProcessNames()
%RETRIEVEPROCESSNAMES This function retrieves all the names for all the
%classes in a folder. 

%Save the Current director and move to the correct directory for the noise
%filters.
previousLocation = pwd;
[localDirectory,~,~] = fileparts(mfilename('fullpath'));
cd(localDirectory);
cd ..\image_handling\noise_filters;
Files = dir('*.m');
num_files = length(Files);

%Create the array to store the name variables. The first value is set to
%nothing.
nameArray = {};
nameArray{1} = 'None';

%Iterate through all the files, ignore the template abstract class. The
%name values are put into a list to be returned.
for i=1:num_files
    [pathstr, name, ext] = fileparts(Files(i).name);
    if strcmp(name,'NoiseFilterTemplate') == false
        nameArray{end+1} = feval(name).processName;
    end
end
cd(previousLocation);
output = nameArray;

end

