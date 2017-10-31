function [ New_matrix ] = Reorder_Matrix(Matrix, map)
    
    New_matrix=zeros(size(Matrix));
    
    nnd=length(map);

    if nnd ~= size(Matrix,1)/3
        error('matrix and map does not match');
    end
    
    for ri=1:nnd
        row_id=map(ri);
        for rd=1:3
           row=(row_id-1)*3+rd;
           
           for ci=1:nnd
               col_id=map(ci);
               for cd=1:3
                   col=(col_id-1)*3+cd;
                   
                   New_matrix((ri-1)*3+rd, (ci-1)*3+cd)=Matrix(row, col);
               end
           end
           
        end
    end

end

