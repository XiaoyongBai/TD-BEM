%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Function: Assemble element matrices

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M_new] = AssembleMK(M_old, M_element, Ele_ids)

    M_new=M_old;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assemble G_element and H_element     
    nnd=length(Ele_ids);

    for a=1:nnd
        
        id_row=Ele_ids(a);
        
        for b=1:3
            local_row=(a-1)*3+b;
            global_row=(id_row-1)*3+b;
            
            for c=1:nnd
               id_col=Ele_ids(c);
                for d=1:3
                    local_col=(c-1)*3+d;
                    global_col=(id_col-1)*3+d;
                    
                    M_new(global_row, global_col)=M_new(global_row, global_col)+M_element(local_row, local_col);               
                end
            end
            
        end
    end

end

