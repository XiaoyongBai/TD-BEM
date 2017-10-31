addpath('./FE_Matrices');
addpath('./BE_Matrices');
addpath('./TracTrans');

%%%%%%%%solving parameters
num_step=2000;
dt=2.5e-3;
beta=0.25; %Newmark beta
gamma=0.5; %Newmark Gamma


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read BEM matrices 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[IEN_BE, nel_BE, Nodes_BE, nnd_BE]  = Geometry_BEM_lin9();
[IEN_FE, nel_FE, Nodes_FE, nnd_FE, load_FE]  = Geometry_FEM_SurfaceBuilding();

[Interface_BE, Interface_FE]=Find_Interface_FE_BE(Nodes_BE, Nodes_FE);

nnd_interface=length(Interface_BE);
interface_DOF_BE=zeros(1,3*nnd_interface);
for n_i=1:nnd_interface
    for d_i=1:3
        interface_DOF_BE((n_i-1)*3+d_i)=(Interface_BE(n_i)-1)*3+d_i;
    end
end
    
[ Matrix_T ] = GlobalT( IEN_BE, Nodes_BE);
T=Matrix_T(interface_DOF_BE, interface_DOF_BE);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read finite element matrices 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load_FE=load_FE*0.1;

%reodering the interfacial nodes of FE zone
[Interface_FE] = Reorder_FE_Interface( Interface_FE, Nodes_FE, Interface_BE, Nodes_BE);

%generate mapping vector
%FE_map(i) is the node which occupy position i in the new matrix
%FE_map_inverse(i) is the new position of node i
[FE_map, FE_map_inverse]= MapVector( nnd_FE, Interface_FE);
[BE_map, BE_map_inverse]= MapVector( nnd_BE, Interface_BE);


%reorder matrices
[M, K] = Read_MK();

[M] = Reorder_Matrix(M, FE_map);
[K] = Reorder_Matrix(K, FE_map);
[load_FE]= Reorder_Vector(load_FE, FE_map);

%form left hand side of FEM
K=K+M/(beta*dt^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_FE_dof=nnd_FE*3;
num_BE_dof=nnd_BE*3;
num_interface_dof=nnd_interface*3;

global_u_dof=num_FE_dof+num_BE_dof-num_interface_dof; %displacement dof numbers
global_t_dof=num_interface_dof; %traction dof numumbers
global_dof=global_u_dof+global_t_dof;

FE_a=zeros(num_step, num_FE_dof);
FE_v=zeros(num_step, num_FE_dof);
FE_u=zeros(num_step, num_FE_dof);

FE_a_original=zeros(num_step, num_FE_dof);
FE_v_original=zeros(num_step, num_FE_dof);
FE_u_original=zeros(num_step, num_FE_dof);

BE_u=zeros(num_step, num_BE_dof);
BE_t=zeros(num_step, num_BE_dof);

%start the velocity at time 0
a0=linsolve(M, load_FE');
FE_a(1,:)=transpose(a0);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num_matrix, H1_DB, H2_DB, G1_DB, G2_DB ] = ReadGH();

H_DB=cell(num_matrix+1,1);
G_DB=cell(num_matrix+1,1);

H_DB_tail=cell(2,1);
G_DB_tail=cell(2,1); 

G_DB{1}=G1_DB{1};
G_DB{2}=G2_DB{1};
H_DB{1}=H1_DB{1};
H_DB{2}=H2_DB{1};

LHS=zeros(global_dof);
RHS=zeros(global_dof,1);

for curr_step=2:num_step    
    
    %form left and right hand side of BEM                
    if curr_step==2
        G_DB{curr_step}=G_DB{curr_step}+G1_DB{curr_step};
        G_DB{curr_step+1}=G2_DB{curr_step};
        H_DB{curr_step}=H_DB{curr_step}+H1_DB{curr_step};
        H_DB{curr_step+1}=H2_DB{curr_step};
        
        G_DB_tail{1}=G2_DB{1};
        G_DB_tail{2}=G2_DB{1};
        H_DB_tail{1}=H2_DB{1};
        H_DB_tail{2}=H2_DB{1};
        
    elseif curr_step>2 && curr_step<=num_matrix
        G_DB{curr_step}=G_DB{curr_step}+G1_DB{curr_step};
        G_DB{curr_step+1}=G2_DB{curr_step};
        H_DB{curr_step}=H_DB{curr_step}+H1_DB{curr_step};
        H_DB{curr_step+1}=H2_DB{curr_step};
        
        G_DB_tail{1}=G2_DB{curr_step-2};
        G_DB_tail{2}=G2_DB{curr_step-1};
        
        H_DB_tail{1}=H2_DB{curr_step-2};
        H_DB_tail{2}=H2_DB{curr_step-1};
    end
    
    [G, H, BE_rhs]=BEM_averaging(G_DB, H_DB, G_DB_tail, H_DB_tail, BE_u, BE_t, curr_step, num_matrix+1);
    
    [G] = Reorder_Matrix(G, BE_map);
    [H] = Reorder_Matrix(H, BE_map);
    [BE_rhs]= Reorder_Vector(BE_rhs, BE_map);
    
    
    %form right hand side of FEM
    u_last=transpose(FE_u(curr_step-1,:));
    v_last=transpose(FE_v(curr_step-1,:));
    a_last=transpose(FE_a(curr_step-1,:));
    
    FE_rhs=transpose(load_FE);
    if (curr_step>50) 
        FE_rhs=FE_rhs*0;
    end
    FE_rhs=FE_rhs+M*(u_last+dt*v_last+0.5*dt^2*(1-2*beta)*a_last)/(dt^2);
    
    
    %Group global lhs and rhs
    inter_FE=num_FE_dof-num_interface_dof;
    inter_BE=num_BE_dof-num_interface_dof;
    
    LHS(1:inter_FE, 1:inter_FE)=K(1:inter_FE, 1:inter_FE);
    LHS(1:inter_FE, inter_FE+inter_BE+1:global_u_dof)=K(1:inter_FE, inter_FE+1:num_FE_dof);
    
    BE_amplifi=1e12;
    LHS(inter_FE+1:global_u_dof, inter_FE+1:global_u_dof)=H(:,:)*BE_amplifi;
    LHS(inter_FE+1:global_u_dof,global_u_dof+1:global_dof)=-G(:,inter_BE+1:num_BE_dof)*BE_amplifi;
    
    LHS(global_u_dof+1:global_dof, 1:inter_FE)=K(inter_FE+1:num_FE_dof, 1:inter_FE);
    LHS(global_u_dof+1:global_dof, inter_FE+inter_BE+1:global_u_dof)=K(inter_FE+1:num_FE_dof,  inter_FE+1:num_FE_dof);
    
    LHS(global_u_dof+1:global_dof,global_u_dof+1:global_dof)=T;
    
    
    RHS(1:inter_FE)=FE_rhs(1:inter_FE);
    RHS(inter_FE+1:global_u_dof)=BE_rhs*BE_amplifi;
    RHS(global_u_dof+1:global_dof)=FE_rhs(inter_FE+1:num_FE_dof);
    
    result=linsolve(LHS, RHS);
    
    FE_u(curr_step, 1:inter_FE)=result(1:inter_FE);
    FE_u(curr_step, inter_FE+1:num_FE_dof)=result(inter_FE+inter_BE+1:global_u_dof);
    
    BE_u_temp=result(inter_FE+1:global_u_dof);
    BE_u(curr_step,:)=Reorder_Vector(BE_u_temp, BE_map_inverse);
    
    BE_t_temp=zeros(1,num_BE_dof);
    BE_t_temp(inter_BE+1:num_BE_dof)=result(global_u_dof+1:global_dof);
    BE_t(curr_step,:)=Reorder_Vector(BE_t_temp, BE_map_inverse);
   
    
    %update
    u_curr=transpose(FE_u(curr_step,:));
    a_curr=(u_curr-u_last-dt*v_last-0.5*dt^2*(1-2*beta)*a_last)/(beta*dt^2);
    v_curr=v_last+dt*((1-gamma)*a_last+gamma*a_curr);
    
    FE_a(curr_step,:)=a_curr;
    FE_v(curr_step,:)=v_curr;
    
    u_curr=Reorder_Vector(u_curr, FE_map_inverse);
    v_curr=Reorder_Vector(v_curr, FE_map_inverse);
    a_curr=Reorder_Vector(a_curr, FE_map_inverse);
    FE_u_original(curr_step,:)=u_curr;
    FE_v_original(curr_step,:)=v_curr;
    FE_a_original(curr_step,:)=a_curr;
end


WriteVTK_FEM_lin(num_step, nnd_FE, Nodes_FE, nel_FE, IEN_FE, FE_u_original, FE_v_original, FE_a_original);
WriteVTK_BEM_lin(num_step, nnd_BE, Nodes_BE, nel_BE, IEN_BE, BE_u);


out_id=1;
filename='BE_dis.txt';
fileID = fopen(filename, 'w');
    
fprintf(fileID,'%s \n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fileID,'%s \n', '% soil structure interaction');
fprintf(fileID,'%s \n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    
fprintf(fileID, '%15s %15s %15s %15s\n', '% t', 'Ux', 'Uy', 'Uz');

    for i=1:num_step
        t=(i-1)*dt;
        fprintf(fileID, '%15e %15e %15e %15e \n', ...
                     t, BE_u(i,out_id*3-2), BE_u(i,out_id*3-1), BE_u(i,out_id*3));
    end
fclose(fileID);


out_id=1;
filename='trac.txt';
fileID = fopen(filename, 'w');
    
fprintf(fileID,'%s \n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fileID,'%s \n', '% soil structure interaction');
fprintf(fileID,'%s \n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    
fprintf(fileID, '%15s %15s %15s %15s\n', '% t', 'Tx', 'Ty', 'Tz');

    for i=1:num_step
        t=(i-1)*dt;
        fprintf(fileID, '%15e %15e %15e %15e \n', ...
                     t, BE_t(i,out_id*3-2), BE_t(i,out_id*3-1), BE_t(i,out_id*3));
    end
fclose(fileID);



out_id=93;
filename='FE_dis.txt';
fileID = fopen(filename, 'w');
    
fprintf(fileID,'%s \n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fileID,'%s \n', '% soil structure interaction');
fprintf(fileID,'%s \n', '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    
fprintf(fileID, '%15s %15s %15s %15s\n', '% t', 'Ux', 'Uy', 'Uz');

    for i=1:num_step
        t=(i-1)*dt;
        fprintf(fileID, '%15e %15e %15e %15e \n', ...
                     t, FE_u_original(i,out_id*3-2), FE_u_original(i,out_id*3-1), FE_u_original(i,out_id*3));
    end
fclose(fileID);


figure;
hold on;
plot((0:num_step-1)*dt, BE_u(:,1),'-r');
plot((0:num_step-1)*dt, FE_u_original(:,93*3-2), '-k');
xlim([0,1]);
figure;
plot((0:num_step-1)*dt, BE_t(:,1*3));
xlim([0,1]);
figure;
plot((0:num_step-1)*dt, FE_u_original(:,107*3-2));
xlim([0,5]);
