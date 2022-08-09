function pf2Dinprocess(fileID)
%% Get a total of element types, and then place the cursor at the beginning
Ele = cell(2,2);
Type = cell(2,1);
nEle = [0;0];
numType = 0;
while ~feof(fileID)
    tline = fgetl(fileID);
    if strncmp(tline,'*Element,',9)
        numType = numType+1;
        %% Read element information
        Type{numType} = tline(16:end);%Read the element type and store it
        
        num_line=ftell(fileID);
        tline = fgetl(fileID);
        tElementNode = textscan(tline,'%f','Delimiter',',');
        fseek(fileID,num_line,'bof');    
        %Adjust the cursor position, ans = 0 means the operation is successful
        EleNode = textscan(fileID,repmat('%f ',1,numel(tElementNode{1})),'Delimiter',',');
        num_line = ftell(fileID);
        
        EleNode = [EleNode{:}]; 
        numEle = size(EleNode,1);
        EleNode2 = EleNode;
        EleNode2(:,1) = EleNode2(:,1)+1e6;
        
        Ele{numType,1} = EleNode;
        Ele{numType,2} = EleNode2;
        nEle(numType) = numEle;
        clear tElementNode tline EleNode EleNode2 numEle
        
    end
end
%Information is written to three log files
%nEle
%numType
%Ele
frewind(fileID);
clear tline

%% Create or overwrite the temp.inp file
%The updated information will be stored here
%The default location is the current folder, which can be changed here
%If it is in the same directory, it must be distinguished from the original file name
fileID2 = fopen('fuck.inp','w+','n','GB2312');

%% Write the beginning of the file
% Read the *element, and write the previous information to temp.inp
tline = fgets(fileID);
fprintf(fileID2,tline);
while 1
    tline = fgets(fileID);
    if strncmp(tline,'*Element,',9)
        break
    end
    fprintf(fileID2,tline);
end

clear tline

%% Write element information

pCPE4 = find(strcmp(Type,'CPE4'));
pCPE3 = find(strcmp(Type,'CPE3'));
% The first layer of CPE4
if ~isempty(pCPE4)
        % Write the first header
        %=================
        fprintf(fileID2,'*************************************************\n');
        fprintf(fileID2,'*User element, nodes=4, type=U1, properties=5, coordinates=2, var=40\n');
        fprintf(fileID2,'1,2\n1,3\n');
        fprintf(fileID2,'*************************************************\n');
        fprintf(fileID2,'*ELEMENT, TYPE=U1, ELSET=SOLID\n');
        % Write data
        fprintf(fileID2,'%d,\t%d,\t%d,\t%d,\t%d\n',Ele{2,1}');
end
        
% The first layer of CPE3
if ~isempty(pCPE3)
        % Write the first header
        %==================
        fprintf(fileID2,'*************************************************\n');
        fprintf(fileID2,'*User element, nodes=3, type=U2, properties=5, coordinates=2, var=30\n');
        fprintf(fileID2,'1,2\n1,3\n');
        fprintf(fileID2,'*************************************************\n');
        fprintf(fileID2,'*ELEMENT, TYPE=U2, ELSET=SOLID\n');
        % Write data
        fprintf(fileID2,'%d,\t%d,\t%d,\t%d\n',Ele{1,1}');
end
% ******************************************
        fprintf(fileID2,'*************************************************\n');
        fprintf(fileID2,'*UEL PROPERTY, ELSET=SOLID\n');
        fprintf(fileID2,'116000., 0.34, 10., 2.7, 1e-07\n');
        fprintf(fileID2,'**E,nu,lc,gc,k\n');
        fprintf(fileID2,'*************************************************\n');
% ******************************************
% The second layer of CPE4
if ~isempty(pCPE4)
        % Write the second header
        %==================
        fprintf(fileID2,'*ELEMENT, TYPE=CPE4, ELSET=Visualization\n');
        % Write data
        fprintf(fileID2,'%d,\t%d,\t%d,\t%d,\t%d\n',Ele{2,2}');
end
        
% The second layer of CPE3
if ~isempty(pCPE3)
        %Write the second header
        %===================
        fprintf(fileID2,'*ELEMENT, TYPE=CPE3, ELSET=Visualization\n');
        % Write data
        fprintf(fileID2,'%d,\t%d,\t%d,\t%d\n',Ele{1,2}');
end
%% Write the remaining data
%And change the corresponding part

fseek(fileID,num_line,'bof');   %Adjust cursor position

while ~feof(fileID)
    tline = fgets(fileID);
    if strncmp(tline,'*Solid Section,',15)
        templ = tline;
        tline = [templ(1:strfind(templ,'elset=')+5),'Visualization',templ(strfind(templ,', material'):end)];
    elseif strncmp(tline,'*User Material,',15)
        fprintf(fileID2,'*Depvar\n');
        fprintf(fileID2,'     9,\n');
        fprintf(fileID2,'1, S11, S11\n');
        fprintf(fileID2,'2, S22, S22\n');
        fprintf(fileID2,'3, S33, S33\n');
        fprintf(fileID2,'4, S12, S12\n');
        fprintf(fileID2,'5, E11, E11\n');
        fprintf(fileID2,'6, E22, E22\n');
        fprintf(fileID2,'7, E33, E33\n');
        fprintf(fileID2,'8, E12, E12\n');
        fprintf(fileID2,'9, PHI, PHI\n');
        tline = '*User Material, constants=0\n';
    end
    fprintf(fileID2,tline);
end

clear tline templ

