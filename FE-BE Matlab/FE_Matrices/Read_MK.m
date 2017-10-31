function [M, K] = Read_MK()

    [ ~, M_E, K_E] = readAnsysEleMCK('MCK_SurfaceBuilding.out');

    [IEN, nel, ~, nnd]  = Geometry_FEM_SurfaceBuilding();
    
    M=zeros(nnd*3);
    K=zeros(nnd*3);
    
    for e_i=1:nel

        M_temp=M_E{e_i};
        K_temp=K_E{e_i};
        
        M=AssembleMK(M, M_temp, IEN(e_i,:));
        K=AssembleMK(K, K_temp, IEN(e_i,:));
        
    end
    
   
    
    

end

