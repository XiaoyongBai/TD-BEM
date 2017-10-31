%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Function: Assemble element matrices

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T_new] = AssembleT(T_old, T_element, Ele_ids)

    T_new=T_old;
    
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
                    
                    T_new(global_row, global_col)=T_new(global_row, global_col)+T_element(local_row, local_col);               
                end
            end
            
        end
    end

end

