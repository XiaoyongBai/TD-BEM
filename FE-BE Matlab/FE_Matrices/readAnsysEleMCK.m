function [ EId, M, K ] = readAnsysEleMCK( file )
    % EId??????? 
    % M ??cell???????? 
    % K ??cell???????? 
    clc; 

    %---------------------------------------- 
    % --- ??????? --- 
    %---------------------------------------- 
    fid = fopen(file, 'r'); % ???? 
    eleCnt = 0; 
    while ~feof(fid) % ?????? 
        tline = fgetl(fid); % ???? 
        if ~isempty(strfind(tline,'STIFFNESS MATRIX FOR ELEMENT')) 
            eleCnt = eleCnt + 1; 
        end 
    end 
    fclose(fid); 


    %---------------------------------------- 
    % --- ???????????????? --- 
    %---------------------------------------- 
    EId = zeros(1,eleCnt); 
    K = cell(1,eleCnt); 
    M = cell(1,eleCnt); 

    fid = fopen(file, 'r'); 

    num = 0; % ??????? 
    n = 0; % ?????? 
    nn = 0; % ???????? 

    while ~feof(fid) 
        tline = fgetl(fid); 

        value = []; 
        if ~isempty(strfind(tline,'STIFFNESS MATRIX FOR ELEMENT')) 
            num = num + 1; 

            idx = regexp(tline,'\d'); 
            eleId = tline(idx); 
            EId(num) = str2num(eleId); 

            tline = fgetl(fid); 

            while ~isempty(tline) 
                value = strcat(value,tline); 
                tline = fgetl(fid); 
            end 

            value = strtrim(value); 
            values = regexp(value, '\s+', 'split'); 

            nn = length(values); 
            for i = 1:nn 
                if isempty(strfind(values{i},'.')) 
                nn = nn - 1; 
                end 
            end 

            n = sqrt(nn); 
            ki = zeros(n,n); 

            cnt = 1; 
            for i = 1:n 
                for j = 1:n 
                    ki(i,j) = str2double(values{cnt+j}); 
                end 
                cnt = cnt + n + 1; 
            end 

            K{num} = ki; 

        elseif ~isempty(strfind(tline,'MASS MATRIX FOR ELEMENT')) 
            value = []; 

            tline = fgetl(fid); 
            while ~isempty(tline) 
                value = strcat(value,tline); 
                tline = fgetl(fid); 
            end 

            value = strtrim(value); 
            values = regexp(value, '\s+', 'split'); 

            mi = zeros(n,n); 

            cnt = 1; 
            for i = 1:n 
                for j = 1:n 
                    mi(i,j) = str2double(values{cnt+j}); 
                end 
            cnt = cnt + n + 1; 
            end 

            M{num} = mi; 
        end 

    end 
    
    fclose(fid); 



end

