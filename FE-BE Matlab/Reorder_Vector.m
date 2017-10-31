function [ New_vector ] = Reorder_Vector(vector, map)
    
    

    New_vector=zeros(1,length(vector));
    
    nnd=length(map);

    if nnd ~= length(vector)/3
        error('vector and map does not match');
    end
    
    for n_i=1:nnd
        id=map(n_i);
        for d_i=1:3
           position=(id-1)*3+d_i;
         
           New_vector((n_i-1)*3+d_i)=vector(position);
           
        end
    end

end

