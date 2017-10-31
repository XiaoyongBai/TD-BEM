function [ new_interface_FE] = Reorder_FE_Interface( interface_FE, Nodes_FE, interface_BE, Nodes_BE)
    

    nnd=length(interface_FE);
    new_interface_FE=zeros(1,nnd);
    i_Nodes_FE=zeros(nnd,3);
    i_Nodes_BE=zeros(nnd,3);
    
    %coordinates of the interfacial nodes
    for i=1:nnd
        fe_id=interface_FE(i);
        be_id=interface_BE(i);
        
        i_Nodes_FE(i,:)=Nodes_FE(fe_id,:);
        i_Nodes_BE(i,:)=Nodes_BE(be_id,:);       
    end
    
    for b_i=1:nnd %loop over boundary element interface
        
        match=0;
        
        for f_i=1:nnd %loop over finite element interface
            dist2=(i_Nodes_BE(b_i,1)-i_Nodes_FE(f_i,1))^2+...
                  (i_Nodes_BE(b_i,2)-i_Nodes_FE(f_i,2))^2+...
                  (i_Nodes_BE(b_i,3)-i_Nodes_FE(f_i,3))^2;
              
            if dist2<1e-10
                match=1;
                new_interface_FE(b_i)=interface_FE(f_i);
                break;
            end
        end
        
        if match==0
            error('FE/BE interface does not match');
        end
        
    end

end

