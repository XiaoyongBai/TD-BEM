%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Function: compute shape value at a specific point
%
%Input: ele_type=type of the element
%       Nodes=coordinates of the element
%       p_1=first shape coordinate
%       p_2=second shape coordinate
%
%Output: N=shape functions
%        X=natural coordinate of the point
%        D_X_1=the first derivative with respect to parametrical coordinate
%        D_X_2=the second derivative with respect to parametrical coordinate
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N, X, D_X_1, D_X_2] = Shape_Value(ele_type, Nodes, p_1, p_2)

    switch ele_type
        case 1
            
            %%%%%%%%%
            % compute shape function   
            N(1)=0.25*(1-p_1)*(1-p_2);
            N(2)=0.25*(1+p_1)*(1-p_2);
            N(3)=0.25*(1+p_1)*(1+p_2);
            N(4)=0.25*(1-p_1)*(1+p_2);

            X=zeros(1,3);
            for a=1:4
                for d=1:3
                    X(d)=X(d)+N(a)*Nodes(a,d);
                end
            end

            %%%%%%%%%
            % compute derivative of shape function
            DN_1st(1)=-0.25*(1-p_2);
            DN_1st(2)=0.25*(1-p_2);
            DN_1st(3)=0.25*(1+p_2);
            DN_1st(4)=-0.25*(1+p_2);


            DN_2nd(1)=-0.25*(1-p_1);
            DN_2nd(2)=-0.25*(1+p_1);
            DN_2nd(3)=0.25*(1+p_1);
            DN_2nd(4)=0.25*(1-p_1);   

            %%%%%%%%%
            % compute directional derivative
            D_X_1=zeros(2,3);

            for a=1:4
                for d=1:3
                    D_X_1(1,d)=D_X_1(1,d)+DN_1st(a)*Nodes(a,d);
                    D_X_1(2,d)=D_X_1(2,d)+DN_2nd(a)*Nodes(a,d);
                end
            end    


            %%%%%%%%
            % compute second order derivative with respect to parametric coordinate
            DDN_1=[0, 0.25, 0.25,0];
            DDN_2=[0,-0.25,-0.25,0];
            DDN_3=[0, 0.25, 0.25,0];
            DDN_4=[0,-0.25,-0.25,0];

            D_X_2=zeros(4,3);

            for d=1:3
                D_X_2(1,d)=D_X_2(1,d)+DDN_1(1)*Nodes(1,d)+DDN_2(1)*Nodes(2,d)+DDN_3(1)*Nodes(3,d)+DDN_4(1)*Nodes(4,d);
                D_X_2(2,d)=D_X_2(2,d)+DDN_1(2)*Nodes(1,d)+DDN_2(2)*Nodes(2,d)+DDN_3(2)*Nodes(3,d)+DDN_4(2)*Nodes(4,d);
                D_X_2(3,d)=D_X_2(3,d)+DDN_1(3)*Nodes(1,d)+DDN_2(3)*Nodes(2,d)+DDN_3(3)*Nodes(3,d)+DDN_4(3)*Nodes(4,d);
                D_X_2(4,d)=D_X_2(4,d)+DDN_1(4)*Nodes(1,d)+DDN_2(4)*Nodes(2,d)+DDN_3(4)*Nodes(3,d)+DDN_4(4)*Nodes(4,d);
            end
            
        case 2
            N(1) = 0.25*(1 - p_1)*(1 - p_2)*(-p_1 - p_2 - 1);
            N(2) = 0.25*(1 + p_1)*(1 - p_2)*(p_1 - p_2 - 1);
            N(3) = 0.25*(1 + p_1)*(1 + p_2)*(p_1 + p_2 - 1);
            N(4) = 0.25*(1 - p_1)*(1 + p_2)*(-p_1 + p_2 - 1);
            N(5) = 0.5*(1 - p_1^2)*(1 - p_2);
            N(6) = 0.5*(1 + p_1)*(1 - p_2^2);
            N(7) = 0.5*(1 - p_1^2)*(1 + p_2);
            N(8) = 0.5*(1 - p_1)*(1 - p_2^2);
            
            X=zeros(1,3);
            for a=1:8
                for d=1:3
                    X(d)=X(d)+N(a)*Nodes(a,d);
                end
            end
            
            DN_1st(1)= -0.25*(1 - p_1)*(1 - p_2) - 0.25*(1 - p_2)*(-1 - p_1 - p_2);
            DN_1st(2)= 0.25*(1 + p_1)*(1 - p_2) + 0.25*(1 - p_2)*(-1 + p_1 - p_2);
            DN_1st(3)=0.25*(1 + p_1)*(1 + p_2) + 0.25*(1 + p_2)*(-1 + p_1 + p_2);
            DN_1st(4)=-0.25*(1 - p_1)*(1 + p_2) - 0.25*(1 + p_2)*(-1 - p_1 + p_2);
            DN_1st(5)=-p_1*(1 - p_2);
            DN_1st(6)=0.5*(1 - p_2^2);
            DN_1st(7)=-p_1*(1 + p_2);
            DN_1st(8)=-0.5*(1 - p_2^2);
            
            DN_2nd(1)=-0.25*(1 - p_1)*(1 - p_2) - 0.25*(1 - p_1)*(-1 - p_1 - p_2);
            DN_2nd(2)=-0.25*(1 + p_1)*(1 - p_2) - 0.25*(1 + p_1)*(-1 + p_1 - p_2);
            DN_2nd(3)=0.25*(1 + p_1)*(1 + p_2) + 0.25*(1 + p_1)*(-1 + p_1 + p_2);
            DN_2nd(4)=0.25*(1 - p_1)*(1 + p_2) + 0.25*(1 - p_1)*(-1 - p_1 + p_2);
            DN_2nd(5)=-0.5*(1 - p_1^2);
            DN_2nd(6)=-(1 + p_1)*p_2;
            DN_2nd(7)=0.5*(1 - p_1^2);
            DN_2nd(8)=-(1 - p_1)*p_2;
            
            %%%%%%%%%
            % compute directional derivative
            D_X_1=zeros(2,3);
            
            for a=1:8
                for d=1:3
                    D_X_1(1,d)=D_X_1(1,d)+DN_1st(a)*Nodes(a,d);
                    D_X_1(2,d)=D_X_1(2,d)+DN_2nd(a)*Nodes(a,d);
                end
            end  
            
            %%%%%%%%
            % compute second order derivative with respect to parametric coordinate
            % DDN_1 contains the four second order derivatives of N1, while
            % DDN_2 contains the four second order derivatives of N2
            DDN_1=[0.5*(1-p_2), 0.25*(1-p_1)+0.25*(1-p_2)+0.25*(-1-p_1-p_2), 0.25*(1-p_1)+0.25*(1-p_2)+0.25*(-1-p_1-p_2), 0.5*(1-p_1)];
            DDN_2=[0.5*(1-p_2), -0.25*(1+p_1)-0.25*(1-p_2)-0.25*(-1+p_1-p_2), -0.25*(1+p_1)-0.25*(1-p_2)-0.25*(-1+p_1-p_2), 0.5*(1+p_1)];
            DDN_3=[0.5*(1+p_2), 0.25*(1+p_1)+0.25*(1+p_2)+0.25*(-1+p_1+p_2), 0.25*(1+p_1)+0.25*(1+p_2)+0.25*(-1+p_1+p_2), 0.5*(1+p_1)];
            DDN_4=[0.5*(1+p_2), -0.25*(1-p_1)-0.25*(1+p_2)-0.25*(-1-p_1+p_2), -0.25*(1-p_1)-0.25*(1+p_2)-0.25*(-1-p_1+p_2), 0.5*(1-p_1)];
            DDN_5=[-(1-p_2), p_1, p_1, 0];
            DDN_6=[0, -p_2, -p_2, -(1+p_1)];
            DDN_7=[-(1+p_2), -p_1, -p_1, 0];
            DDN_8=[0, p_2, p_2, -(1-p_1)];

            D_X_2=zeros(4,3);

            for d=1:3
                D_X_2(1,d)=DDN_1(1)*Nodes(1,d)+DDN_2(1)*Nodes(2,d)+DDN_3(1)*Nodes(3,d)+DDN_4(1)*Nodes(4,d)+DDN_5(1)*Nodes(5,d)+DDN_6(1)*Nodes(6,d)+DDN_7(1)*Nodes(7,d)+DDN_8(1)*Nodes(8,d);
                D_X_2(2,d)=DDN_1(2)*Nodes(1,d)+DDN_2(2)*Nodes(2,d)+DDN_3(2)*Nodes(3,d)+DDN_4(2)*Nodes(4,d)+DDN_5(2)*Nodes(5,d)+DDN_6(2)*Nodes(6,d)+DDN_7(2)*Nodes(7,d)+DDN_8(2)*Nodes(8,d);
                D_X_2(3,d)=DDN_1(3)*Nodes(1,d)+DDN_2(3)*Nodes(2,d)+DDN_3(3)*Nodes(3,d)+DDN_4(3)*Nodes(4,d)+DDN_5(3)*Nodes(5,d)+DDN_6(3)*Nodes(6,d)+DDN_7(3)*Nodes(7,d)+DDN_8(3)*Nodes(8,d);
                D_X_2(4,d)=DDN_1(4)*Nodes(1,d)+DDN_2(4)*Nodes(2,d)+DDN_3(4)*Nodes(3,d)+DDN_4(4)*Nodes(4,d)+DDN_5(4)*Nodes(5,d)+DDN_6(4)*Nodes(6,d)+DDN_7(4)*Nodes(7,d)+DDN_8(4)*Nodes(8,d);
            end
            
        otherwise
            error('Shape_value:unsupported element type');
    end

end

