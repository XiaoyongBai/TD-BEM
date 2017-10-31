function [num_matrix, H1_DB, H2_DB, G1_DB, G2_DB ] = ReadGH()

    num_matrix=8;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %read matrices
    H1_DB=cell(num_matrix,1);
    H2_DB=cell(num_matrix,1);
    G1_DB=cell(num_matrix,1);
    G2_DB=cell(num_matrix,1);

    for i=1:num_matrix
        H1_name=sprintf('%s%d%s','./BE_Matrices/SurfaceBuilding_DenseSoil/H1_step',i-1,'.txt');
        H1_DB{i}=load(H1_name);

        H2_name=sprintf('%s%d%s','./BE_Matrices/SurfaceBuilding_DenseSoil/H2_step',i-1,'.txt');
        H2_DB{i}=load(H2_name);

        G1_name=sprintf('%s%d%s','./BE_Matrices/SurfaceBuilding_DenseSoil/G1_step',i-1,'.txt');
        G1_DB{i}=load(G1_name);

        G2_name=sprintf('%s%d%s','./BE_Matrices/SurfaceBuilding_DenseSoil/G2_step',i-1,'.txt');
        G2_DB{i}=load(G2_name);
    end


end

