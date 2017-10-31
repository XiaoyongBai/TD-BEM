%compute T matrix for boundary element
%
%Nodes=coordinates of the element nodes
%
function [ elment_T ] = ElementT( Nodes )

    %!!this element is only for 4-node linear element
    nnd=4;
    elment_T=zeros(4*3);
    
    nq_1=6;
    nq_2=6;
    
    for q_1=1:nq_1
        for q_2=1:nq_2
            [~, ~, N, w, J] = Shape_Function(1, Nodes, nq_1, nq_2, q_1, q_2);
            
            Shape_matrix=[];
            for node=1:nnd
                matrix_temp=[N(node),0,0;0,N(node),0;0,0,N(node)];
                Shape_matrix=[Shape_matrix,matrix_temp];
            end
            
            elment_T = elment_T + transpose(Shape_matrix)*Shape_matrix*w*J;
        end
    end


end

