function [IEN, nel, Nodes, nnd, interface]  = Geometry_BEM( )
    

    nel=1;
    nnd=4;
    Nodes=[  0.00000   0.00000   0.0
   1.0	     0.0       0.0
   1.0       1.0       0.0
   0.0       1.0       0.0
   
      ];
      
        
   IEN=[        
     1,2,3,4
    ];


    interface=[1,2,3,4];
            


end

