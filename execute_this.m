%% Open the original InP file
% =========Change the path here
filepath = '.\';
% =========Change the file name here
filename = '1.inp';

% =================================
fullname = [filepath filename];
[fileID,errmsg] = fopen(fullname,'r','n','GB2312');
while fileID < 0
    disp(errmsg);
    fullname = input('Open file: ', 's');
    [fileID,errmsg] = fopen(fullname,'r');
end

clear errmsg fullname filepath
%% Get the dimension and place the cursor at the beginning
tline = fgetl(fileID);
while ~strncmp(tline,'*Node',5)
    tline = fgetl(fileID);
end
tline = fgetl(fileID);
tNodeCoor = textscan(tline,'%f','Delimiter',',');
dimension = length(tNodeCoor{1})-1;

frewind(fileID);
clear tNodeCoor tline
%% Processing InP files

switch dimension
    case 2
        pf2Dinprocess(fileID)
    case 3
        pf3Dinprocess(fileID)
    otherwise
        warning('Dimension is incorrect!')
end

fclose all;