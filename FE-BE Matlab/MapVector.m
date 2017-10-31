function [ map, map_inverse ] = MapVector( nnd, interface )
    
    map=zeros(1,nnd);
    map_inverse=zeros(1,nnd);
    
    i_nnd=length(interface); %number of nodes on interface
    
    count=0;
    
    for i=1:nnd
        
        is_interface=0;
        
        for j=1:i_nnd
            if i==interface(j)
                is_interface=1;
            end
        end
        
        if is_interface==0
            count=count+1;
            map(count)=i;
            
            map_inverse(i)=count;
        end
    end

    map(count+1:nnd)=interface;
    map_inverse(interface)=count+1:nnd;
    
    

end

