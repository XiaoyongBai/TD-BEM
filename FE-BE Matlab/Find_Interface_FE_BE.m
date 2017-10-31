%%determine the interfacial nodes between BEM and FEM
%
%input: 
%   BE_Nodes=coordinates of BE nodes
%   FE_Nodes=coordinates of FE nodes
%output:
%   Interface_BE_node=id of 
function [Interface_BE_node, Interface_FE_node]=Find_Interface_FE_BE(BE_Nodes, FE_Nodes)

    Interface_BE_node=[];
    Interface_FE_node=[];
    
    FE_nnd=size(FE_Nodes, 1);
    BE_nnd=size(BE_Nodes, 1);
    
    FE_flag=zeros(FE_nnd, 1);
    BE_flag=zeros(BE_nnd, 1);
        
    for i=1:FE_nnd
        coord_fe=FE_Nodes(i,:);
        
        for j=1:BE_nnd
            coord_be=BE_Nodes(j,:);
            
            r=sqrt( (coord_fe(1)-coord_be(1))^2+(coord_fe(2)-coord_be(2))^2+(coord_fe(3)-coord_be(3))^2 );
            
            if (r<1e-8) %the two points have the same coordinate
               
                if (FE_flag(i)==0) %To avoid multiple add of a FE node
                    Interface_FE_node=[Interface_FE_node, i];
                end
                
                if (BE_flag(j)==0) %To avoid multiple add of a BE node
                    Interface_BE_node=[Interface_BE_node, j];
                end
                
                FE_flag(i)=1;
                BE_flag(j)=1;
            end
            
        end
        
    end

end

