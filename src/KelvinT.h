//
//  KelvinT.h
//  GreenFunction
//
//  Created by Xiaoyong Bai on 5/13/16.
//
//

#ifndef KelvinT_hpp
#define KelvinT_hpp

#include <iostream>

namespace GreenFunction{
    
class KelvinT
{
public:
    
    KelvinT();
    ~KelvinT();
        
    //set material constants
    //
    //input:
    //		E=Young's modulus
    //      nu=Possion's ratio
    // 		density=density
    void SetMaterial(double E, double nu);
    
    //set positions of source and observation
    //
    //input:
    //      src=coordinate of source point
    //      dest=coordinate of observation point
    //      normal=normal vector of the observation surface
    void SetGeometry(double* src, double* dest, double* normal);

    //drive the computation after all the settings
    void Compute();
    void ComputeDisplacement();
    void ComputeStress();
    void ComputeTraction();
    
    
    
    //return 1 if i==j, otherwise return 0
    //
    int KronDelta(int i, int j);
    
    //Get pointer to the Green's functions.
    //For displacement and tractions, the first index stands for component, and 2nd is loading direction.
    //SX is stress corresponding to loading in X-direction. Similar for Y and Z.
    //
    const double* GetTraction();
    const double* GetDisplacement();
    const double* GetStress(const double * SX, const double * SY, const double * SZ);

    
        
private:

    // fE=Young's modulus
    // fNu=Possion's ratio
    // fMu=Shear modulus
    // fDensity=density
    // fR_v=directional vector from source point to observation point
    // fR=distance from source point to observation point
    // fN=normal of the surface on which the traction is computing
    double fE, fNu, fMu;
    double *fR_v;
    double fR;
    double *fN;
        
        
    //fDis=Kelvin displacement
    //fTrac=Kelvin traction
    //fStress_*=Kelvin stress due to loading in direction *
    double* fDis;
    double* fTrac;
    double* fStress_X;
    double* fStress_Y;
    double* fStress_Z;
    
};//end of definition of class KelvinT
    


    
}//end of namespace



#endif /* KelvinT_hpp */
