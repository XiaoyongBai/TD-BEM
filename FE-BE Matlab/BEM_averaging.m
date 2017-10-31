%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Function: Solve the static problem after obtainning G and H
%
%Input: Problem=problem set
%       C_vector=vector of free term coefficient
%       G_DB=Data base for G matrices
%       H_DB=Data base for H matrices
%       Dis_DB=Database for displacement
%       Trac_Db=Databse for traction
%       curr_step=current step number
%
%Output: Dis_curr=solution of displacement
%        Trac_curr=solution of tractions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G_curr, H_curr, rhs]=BEM_averaging(G_DB, H_DB,  G_DB_tail, H_DB_tail, Dis_DB, Trac_DB, curr_step, max_step)

    if curr_step==2
        G_curr=G_DB{1};
        H_curr=H_DB{1};
        
        rhs=zeros(size(G_curr, 1), 1);   
                
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% specify the main matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        G_curr=4*G_DB{1}+G_DB{2};
        H_curr=4*H_DB{1}+H_DB{2};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% compute rhs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        rhs=zeros(size(G_curr, 1), 1);
        
        %%%add terms by loop
        %n+1 step
        for f_i=3:1:min(max_step,curr_step+1)
            Trac_temp=transpose(Trac_DB(curr_step+2-f_i,:));
            Dis_temp=transpose(Dis_DB(curr_step+2-f_i,:));

            rhs=rhs+G_DB{f_i}*Trac_temp;
            rhs=rhs-H_DB{f_i}*Dis_temp;
        end   
        
        %n step
        for f_i=2:1:min(max_step,curr_step-1)
            Trac_temp=transpose(Trac_DB(curr_step+1-f_i,:));
            Dis_temp=transpose(Dis_DB(curr_step+1-f_i,:));

            rhs=rhs+2*G_DB{f_i}*Trac_temp;
            rhs=rhs-2*H_DB{f_i}*Dis_temp;
        end 
        if curr_step <= max_step
            Trac_temp=transpose(Trac_DB(1,:));
            Dis_temp=transpose(Dis_DB(1,:));

            rhs=rhs+2*G_DB_tail{2}*Trac_temp;
            rhs=rhs-2*H_DB_tail{2}*Dis_temp;
        end
        
        %n-1 step
        for f_i=2:1:min(max_step,curr_step-2)
            Trac_temp=transpose(Trac_DB(curr_step-f_i,:));
            Dis_temp=transpose(Dis_DB(curr_step-f_i,:));

            rhs=rhs+G_DB{f_i}*Trac_temp;
            rhs=rhs-H_DB{f_i}*Dis_temp;
        end
        
        if curr_step-1 <= max_step
            Trac_temp=transpose(Trac_DB(1,:));
            Dis_temp=transpose(Dis_DB(1,:));

            rhs=rhs+G_DB_tail{1}*Trac_temp;
            rhs=rhs-H_DB_tail{1}*Dis_temp;
        end
        
              
    end   
    
    


end

