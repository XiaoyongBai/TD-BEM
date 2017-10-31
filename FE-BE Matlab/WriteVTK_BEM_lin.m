function [ ] = WriteVTK_BEM_lin(num_step, nnd, nodes, nel, IEN, dis_DB)

    for si=1:num_step

        file_name=sprintf('%s%d%s','./BE_VTK_files/step_',si,'.vtk');
        fileID = fopen(file_name, 'w');

        fprintf(fileID, '%s\n', '# vtk DataFile Version 3.0');
        fprintf(fileID, '%s\n', 'Boundary element result');
        fprintf(fileID, '%s\n', 'ASCII');
        fprintf(fileID, '%s\n', 'DATASET UNSTRUCTURED_GRID');
        fprintf(fileID, '%s  %d  %s\n', 'POINTS', nnd, 'float');

        for ni=1:nnd
            fprintf(fileID, '%15e %15e %15e \n', nodes(ni,1), nodes(ni,2), nodes(ni,3));
        end
        
        fprintf(fileID, '%s %d %d \n', 'CELLS ', nel, nel*5);
        
        for ei=1:nel
            fprintf(fileID, '%2d %8d %8d %8d %8d\n', 4,...
                    IEN(ei,1)-1, IEN(ei,2)-1, IEN(ei,3)-1, IEN(ei,4)-1);
        end
        
        fprintf(fileID, '%s %d\n', 'CELL_TYPES ', nel);
        
        for ei=1:nel
            fprintf(fileID, '%d\n', 9);
        end
              
        
        fprintf(fileID, '%s %d\n', 'POINT_DATA ', nnd);
		fprintf(fileID, '%s %s %s\n', 'VECTORS', 'Displacement', 'double');
        dis=dis_DB(si,:);
        for ni=1:nnd
            fprintf(fileID, '%15e %15e %15e \n', dis((ni-1)*3+1), dis((ni-1)*3+2), dis((ni-1)*3+3));
        end


        fclose(fileID);

    end
        



end

