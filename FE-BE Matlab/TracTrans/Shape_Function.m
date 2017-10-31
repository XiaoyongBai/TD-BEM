%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Function: shape function 
%
%Input: ele_type=type of the element
%       Nodes=node coordinates of the element
%       nq_1=number of quadrature points in direction 1
%       nq_2=number of quadrature points in direction 2
%       q_1=current quadrature point in direction 1
%       q_2=current quadrature point in direction 2
%
%Output: 
%        x=coordinates of current quadrature point 
%        n=out normal at current quadrature point
%        N=shape function values
%        w=weight of current quadrature point
%        J=determinant of Jacobian matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, n, N, w, J] = Shape_Function(ele_type, Nodes, nq_1, nq_2, q_1, q_2)

    if (nq_1 ~= nq_2)
        error('Shape_Function: nq_1 != nq_2');
    end
    
    switch nq_1
        case 6
            shape_coords=[  0.6612093864662645;
                            -0.6612093864662645;
                            -0.2386191860831969;
                            0.2386191860831969;
                            -0.9324695142031521;
                            0.9324695142031521];
                        
            weights=[   0.3607615730481386;	
                        0.3607615730481386;	
                        0.4679139345726910;	
                        0.4679139345726910;	
                        0.1713244923791704;	
                        0.1713244923791704];	


        case 10
            shape_coords= [-0.148874339;
                            0.148874339;
                            -0.433395394;
                            0.433395394;
                            -0.679409568;
                            0.679409568;
                            -0.865063367;
                            0.865063367;
                            -0.973906529;
                            0.973906529];
             weights=[0.295524225
                        0.295524225
                        0.269266719
                        0.269266719
                        0.219086363
                        0.219086363
                        0.149451349
                        0.149451349
                        0.066671344
                        0.066671344
                        ];
           

        case 20
            shape_coords=[ -0.076526521;
                            0.076526521;
                            -0.227785851;
                            0.227785851;
                            -0.373706089;
                            0.373706089;
                            -0.510867002;
                            0.510867002;
                            -0.636053681;
                            0.636053681;
                            -0.746331906;
                            0.746331906;
                            -0.839116972;
                            0.839116972;
                            -0.912234428;
                            0.912234428;
                            -0.963971927;
                            0.963971927;
                            -0.993128599;
                            0.993128599
                            ];
            weights=[   0.152753387;
                        0.152753387;
                        0.149172986;
                        0.149172986;
                        0.142096109;
                        0.142096109;
                        0.131688638;
                        0.131688638;
                        0.118194532;
                        0.118194532;
                        0.10193012;
                        0.10193012;
                        0.083276742;
                        0.083276742;
                        0.062672048;
                        0.062672048;
                        0.04060143;
                        0.04060143;
                        0.017614007;
                        0.017614007
                        ];
        otherwise
            error('Shape_Function:unsupported quadrature number');
    end

    %%%%%%%%%
    % compute shape function
    p_1=shape_coords(q_1);
    p_2=shape_coords(q_2);
    
    switch ele_type
        case 1
            N(1)=0.25*(1-p_1)*(1-p_2);
            N(2)=0.25*(1+p_1)*(1-p_2);
            N(3)=0.25*(1+p_1)*(1+p_2);
            N(4)=0.25*(1-p_1)*(1+p_2);
            
            x=zeros(1,3);
            for a=1:4
                for d=1:3
                    x(d)=x(d)+N(a)*Nodes(a,d);
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
            DX_1st=zeros(1,3);
            DX_2nd=zeros(1,3);

            for a=1:4
                for d=1:3
                    DX_1st(d)=DX_1st(d)+DN_1st(a)*Nodes(a,d);
                    DX_2nd(d)=DX_2nd(d)+DN_2nd(a)*Nodes(a,d);

                end
            end    
    
            %%%%%%%%
            % compute outer normal of the surface element
            n=cross(DX_1st, DX_2nd);
            J=norm(n);
            n=n/J;          
            
        case 2
            N(1) = 0.25*(1 - p_1)*(1 - p_2)*(-p_1 - p_2 - 1);
            N(2) = 0.25*(1 + p_1)*(1 - p_2)*(p_1 - p_2 - 1);
            N(3) = 0.25*(1 + p_1)*(1 + p_2)*(p_1 + p_2 - 1);
            N(4) = 0.25*(1 - p_1)*(1 + p_2)*(-p_1 + p_2 - 1);
            N(5) = 0.5*(1 - p_1^2)*(1 - p_2);
            N(6) = 0.5*(1 + p_1)*(1 - p_2^2);
            N(7) = 0.5*(1 - p_1^2)*(1 + p_2);
            N(8) = 0.5*(1 - p_1)*(1 - p_2^2);
            
            x=zeros(1,3);
            for a=1:8
                for d=1:3
                    x(d)=x(d)+N(a)*Nodes(a,d);
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
            DX_1st=zeros(1,3);
            DX_2nd=zeros(1,3);

            for a=1:8
                for d=1:3
                    DX_1st(d)=DX_1st(d)+DN_1st(a)*Nodes(a,d);
                    DX_2nd(d)=DX_2nd(d)+DN_2nd(a)*Nodes(a,d);

                end
            end    
    
            %%%%%%%%
            % compute outer normal of the surface element
            n=cross(DX_1st, DX_2nd);
            J=norm(n);
            n=n/J; 
            
        otherwise
            error('Shape_function:unsupported element type');
    end


    
   

    w=weights(q_1)*weights(q_2);

end

