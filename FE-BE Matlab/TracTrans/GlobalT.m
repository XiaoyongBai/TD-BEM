%compute traction-nodal force transformation matrix for the boundary
function [ Matrix_T ] = GlobalT( IEN, Nodes )

    nel=size(IEN,1);
    nnd_ele=size(IEN,2);
    nnd=size(Nodes,1);
    
    Matrix_T=zeros(nnd*3);
    
    for e_i=1:nel
        ele_nodes=zeros(nnd_ele, 3);
    
        for a=1:nnd_ele
            ele_nodes(a, 1:3)=Nodes(IEN(e_i, a), 1:3);
        end
        
        element_T=ElementT(ele_nodes);
        
        Matrix_T=AssembleT(Matrix_T, element_T, IEN(e_i,:) );
    end



end

