function pf3Dinprocess(fileID)
%% Get a total of element types, and then place the cursor at the beginning
Ele = cell(4,2);
Type = cell(4,1);
nEle = [0;0;0;0];
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
fileID2 = fopen('temp3D.inp','w+','n','GB2312');
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

pC3D4 = find(strncmp(Type,'C3D4',4));
pC3D6 = find(strncmp(Type,'C3D6',4));
pC3D8 = find(strncmp(Type,'C3D8',4));
pC3D10 = find(strncmp(Type,'C3D10',4));
% =================================================
% The first layer of C3D4
if ~isempty(pC3D4)
        % Write the first header
        %=======================
        fprintf(fileID2,'*************************************************\n');
        fprintf(fileID2,'*User element, nodes=4, type=U1, properties=5, coordinates=3, var=56\n');
        fprintf(fileID2,'1,2,3\n1,4\n');
        fprintf(fileID2,'*************************************************\n');
        fprintf(fileID2,'*ELEMENT, TYPE=U1, ELSET=SOLID\n');
        % Write data
        fprintf(fileID2,'\t%d,\t%d,\t%d,\t%d,\t%d\n',Ele{pC3D4,1}');
end
        
% The first layer of C3D6
if ~isempty(pC3D6)
        % Write the first header
        %=======================
        fprintf(fileID2,'*************************************************\n');
        fprintf(fileID2,'*User element, nodes=6, type=U2, properties=5, coordinates=3, var=84\n');
        fprintf(fileID2,'1,2,3\n1,4\n');
        fprintf(fileID2,'*************************************************\n');
        fprintf(fileID2,'*ELEMENT, TYPE=U2, ELSET=SOLID\n');
        % Write data
        fprintf(fileID2,'\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d\n',Ele{pC3D6,1}');
end

% The first layer of C3D8
if ~isempty(pC3D8)
        % Write the first header
        %=======================
        fprintf(fileID2,'*************************************************\n');
        fprintf(fileID2,'*User element, nodes=8, type=U3, properties=5, coordinates=3, var=112\n');
        fprintf(fileID2,'1,2,3\n1,4\n');
        fprintf(fileID2,'*************************************************\n');
        fprintf(fileID2,'*ELEMENT, TYPE=U3, ELSET=SOLID\n');
        % Write data
        fprintf(fileID2,'\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d\n',Ele{pC3D8,1}');
end

% The first layer of C3D10
if ~isempty(pC3D10)
        % Write the first header
        %=======================
        fprintf(fileID2,'*************************************************\n');
        fprintf(fileID2,'*User element, nodes=10, type=U4, properties=5, coordinates=3, var=140\n');
        fprintf(fileID2,'1,2,3\n1,4\n');
        fprintf(fileID2,'*************************************************\n');
        fprintf(fileID2,'*ELEMENT, TYPE=U4, ELSET=SOLID\n');
        % Write data
        fprintf(fileID2,'\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d\n',Ele{pC3D10,1}');
end
% ******************************************
        fprintf(fileID2,'*************************************************\n');
        fprintf(fileID2,'*UEL PROPERTY, ELSET=SOLID\n');
        fprintf(fileID2,'75000., 0.3, 1., 200., 1e-07\n');
        fprintf(fileID2,'**E,nu,lc,gc,k\n');
        fprintf(fileID2,'*************************************************\n');
% ******************************************
% ===========================================
% ******************************************
% The second layer of C3D4
if ~isempty(pC3D4)
        %Write the second header
        %=======================
        fprintf(fileID2,'*ELEMENT, TYPE=C3D4, ELSET=Visualization\n');
        % Write data
        fprintf(fileID2,'\t%d,\t%d,\t%d,\t%d,\t%d\n',Ele{pC3D4,2}');
end
 
% The second layer of C3D6
        %Write the second header
        %=======================
        fprintf(fileID2,'*ELEMENT, TYPE=C3D6, ELSET=Visualization\n');
        % Write data
        fprintf(fileID2,'\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d\n',Ele{pC3D6,2}');
end

% The second layer of C3D8
if ~isempty(pC3D8)
        %Write the second header
        %=======================
        fprintf(fileID2,'*ELEMENT, TYPE=C3D8, ELSET=Visualization\n');
        % Write data
        fprintf(fileID2,'\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d\n',Ele{pC3D8,2}');
end

% The second layer of C3D10
if ~isempty(pC3D10)
        %Write the second header
        %=======================
        fprintf(fileID2,'*ELEMENT, TYPE=C3D10, ELSET=Visualization\n');
        % Write data
        fprintf(fileID2,'\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d,\t%d\n',Ele{pC3D10,2}');
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
        fprintf(fileID2,'     13,\n');
        fprintf(fileID2,'1, S11, S11\n');
        fprintf(fileID2,'2, S22, S22\n');
        fprintf(fileID2,'3, S33, S33\n');
        fprintf(fileID2,'4, S12, S12\n');
        fprintf(fileID2,'5, S13, S13\n');
        fprintf(fileID2,'6, S23, S23\n');
        fprintf(fileID2,'7, E11, E11\n');
        fprintf(fileID2,'8, E22, E22\n');
        fprintf(fileID2,'9, E33, E33\n');
        fprintf(fileID2,'10, E12, E12\n');
        fprintf(fileID2,'11, E13, E13\n');
        fprintf(fileID2,'12, E23, E23\n');
        fprintf(fileID2,'13, PHI, PHI\n');
        tline = '*User Material, constants=1\n';
    end
    fprintf(fileID2,tline);
end

clear tline templ
